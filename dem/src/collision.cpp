// Collision detection: spatial hash grid broadphase + narrowphase

#include "collision.hpp"
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dem {

void SpatialGrid::build(const std::vector<Particle>& particles,
                        const ParticleSystem& psys,
                        double cell_size) {
    cell_size_ = cell_size;
    cells_.clear();

    for (uint32_t i = 0; i < (uint32_t)particles.size(); i++) {
        if (!particles[i].active) continue;
        CellKey key = pos_to_cell(particles[i].pos);
        cells_[key].push_back(i);
    }
}

std::vector<CollisionPair> SpatialGrid::find_pairs() const {
    std::vector<CollisionPair> pairs;

    for (auto& [cell_key, indices] : cells_) {
        // Self-cell pairs
        for (size_t i = 0; i < indices.size(); i++) {
            for (size_t j = i + 1; j < indices.size(); j++) {
                pairs.push_back({indices[i], indices[j]});
            }
        }

        // Neighbor cells (only forward neighbors to avoid duplicates)
        for (int dx = 0; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    if (dx == 0 && dy < 0) continue;
                    if (dx == 0 && dy == 0 && dz < 0) continue;

                    CellKey neighbor{cell_key.x + dx, cell_key.y + dy, cell_key.z + dz};
                    auto it = cells_.find(neighbor);
                    if (it == cells_.end()) continue;

                    for (uint32_t a : indices) {
                        for (uint32_t b : it->second) {
                            pairs.push_back({a, b});
                        }
                    }
                }
            }
        }
    }

    return pairs;
}

void CollisionDetector::detect(ParticleSystem& psys) {
    sphere_contacts_.clear();
    wall_contacts_.clear();

    auto& particles = psys.particles;
    if (particles.empty()) return;

    // Determine cell size from max particle radius
    double max_r = 0;
    for (auto& p : particles) {
        if (!p.active) continue;
        double r = psys.bounding_radius(p);
        if (r > max_r) max_r = r;
    }
    double cell_size = 2.1 * max_r;
    if (cell_size < 1e-6) cell_size = 0.01;

    // Broadphase
    grid_.build(particles, psys, cell_size);
    auto pairs = grid_.find_pairs();

    // Narrowphase: sphere-sphere
    for (auto& pair : pairs) {
        auto& a = particles[pair.idx_a];
        auto& b = particles[pair.idx_b];
        if (!a.active || !b.active) continue;

        double ra = psys.bounding_radius(a);
        double rb = psys.bounding_radius(b);

        vec3 diff = b.pos - a.pos;
        double dist = diff.length();
        double overlap = (ra + rb) - dist;

        if (overlap > 0 && dist > 1e-30) {
            SphereContact sc;
            sc.idx_a = pair.idx_a;
            sc.idx_b = pair.idx_b;
            sc.normal = diff / dist;
            sc.overlap = overlap;
            sc.contact_point = a.pos + sc.normal * (ra - overlap * 0.5);
            sphere_contacts_.push_back(sc);
        }
    }

    // Narrowphase: sphere-wall (geometry queries)
    if (geometry_) {
        // Can be parallelized
        #ifdef _OPENMP
        #pragma omp parallel
        {
            std::vector<WallContact> local_contacts;
            #pragma omp for nowait schedule(dynamic, 64)
            for (int i = 0; i < (int)particles.size(); i++) {
                if (!particles[i].active) continue;
                double r = psys.bounding_radius(particles[i]);
                auto result = geometry_->query_sphere(particles[i].pos, r);
                if (result.hit && result.penetration > 0) {
                    WallContact wc;
                    wc.idx = (uint32_t)i;
                    wc.contact_point = result.contact_point;
                    wc.normal = result.contact_normal;
                    wc.overlap = result.penetration;
                    local_contacts.push_back(wc);
                }
            }
            #pragma omp critical
            {
                wall_contacts_.insert(wall_contacts_.end(),
                    local_contacts.begin(), local_contacts.end());
            }
        }
        #else
        for (int i = 0; i < (int)particles.size(); i++) {
            if (!particles[i].active) continue;
            double r = psys.bounding_radius(particles[i]);
            auto result = geometry_->query_sphere(particles[i].pos, r);
            if (result.hit && result.penetration > 0) {
                WallContact wc;
                wc.idx = (uint32_t)i;
                wc.contact_point = result.contact_point;
                wc.normal = result.contact_normal;
                wc.overlap = result.penetration;
                wall_contacts_.push_back(wc);
            }
        }
        #endif
    }
}

} // namespace dem
