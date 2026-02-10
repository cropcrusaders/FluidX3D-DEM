#pragma once
// Collision detection: broadphase (uniform grid) + narrowphase

#include "dem_math.hpp"
#include "particle.hpp"
#include "geometry.hpp"
#include <vector>
#include <unordered_map>

namespace dem {

// Collision pair
struct CollisionPair {
    uint32_t idx_a, idx_b;  // particle indices
};

// Spatial hash grid for broadphase
class SpatialGrid {
public:
    void build(const std::vector<Particle>& particles,
               const ParticleSystem& psys,
               double cell_size);

    // Find all potential sphere-sphere collision pairs
    std::vector<CollisionPair> find_pairs() const;

    double get_cell_size() const { return cell_size_; }

private:
    struct CellKey {
        int32_t x, y, z;
        bool operator==(const CellKey& o) const { return x==o.x && y==o.y && z==o.z; }
    };
    struct CellKeyHash {
        size_t operator()(const CellKey& k) const {
            size_t h = 73856093ULL * (uint32_t)k.x
                     ^ 19349663ULL * (uint32_t)k.y
                     ^ 83492791ULL * (uint32_t)k.z;
            return h;
        }
    };

    double cell_size_ = 0.01;
    std::unordered_map<CellKey, std::vector<uint32_t>, CellKeyHash> cells_;

    CellKey pos_to_cell(const vec3& p) const {
        return {
            (int32_t)std::floor(p.x / cell_size_),
            (int32_t)std::floor(p.y / cell_size_),
            (int32_t)std::floor(p.z / cell_size_)
        };
    }
};

// Full collision detection system
class CollisionDetector {
public:
    void set_geometry(const Geometry* geo) { geometry_ = geo; }

    // Detect all collisions (populates results)
    void detect(ParticleSystem& psys);

    // Results
    struct SphereContact {
        uint32_t idx_a, idx_b;
        vec3 contact_point;
        vec3 normal;  // from a to b
        double overlap;
    };
    struct WallContact {
        uint32_t idx;
        vec3 contact_point;
        vec3 normal;  // outward from wall
        double overlap;
    };

    const std::vector<SphereContact>& sphere_contacts() const { return sphere_contacts_; }
    const std::vector<WallContact>& wall_contacts() const { return wall_contacts_; }

private:
    const Geometry* geometry_ = nullptr;
    SpatialGrid grid_;
    std::vector<SphereContact> sphere_contacts_;
    std::vector<WallContact> wall_contacts_;
};

} // namespace dem
