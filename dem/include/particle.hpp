#pragma once
// Particle / Seed data structures for DEM simulation

#include "dem_math.hpp"
#include <vector>
#include <string>

namespace dem {

enum class SeedShape { SPHERE, MULTI_SPHERE_ELLIPSOID };

// Material properties for contact model
struct Material {
    std::string name = "default";
    double density = 1200.0;       // kg/m^3 (typical seed)
    double youngs_modulus = 1e7;   // Pa (not used for linear spring, but for reference)
    double poisson_ratio = 0.3;
};

// Contact pair properties
struct ContactProperties {
    double kn = 1e4;      // normal stiffness (N/m)
    double kt = 2500.0;   // tangential stiffness (N/m)
    double restitution = 0.5;  // coefficient of restitution
    double friction_s = 0.5;   // static friction coefficient
    double friction_r = 0.01;  // rolling friction coefficient
};

// Sphere component for multi-sphere representation
struct SubSphere {
    vec3 local_pos;  // position in body frame
    double radius;
};

// Seed type definition
struct SeedType {
    std::string name = "cotton";
    SeedShape shape = SeedShape::SPHERE;
    double radius = 0.004;    // m (4mm default)
    double mass = 0.001;      // kg
    // Ellipsoid semi-axes (only used for MULTI_SPHERE_ELLIPSOID)
    double a = 0.004, b = 0.003, c = 0.002;  // semi-axes in m

    // Multi-sphere decomposition (generated from ellipsoid params)
    std::vector<SubSphere> sub_spheres;

    // Inertia tensor (diagonal, body frame)
    vec3 inertia() const {
        if (shape == SeedShape::SPHERE) {
            double I = 0.4 * mass * radius * radius;
            return {I, I, I};
        } else {
            // Ellipsoid inertia
            return {
                0.2 * mass * (b*b + c*c),
                0.2 * mass * (a*a + c*c),
                0.2 * mass * (a*a + b*b)
            };
        }
    }

    // Generate multi-sphere approximation of ellipsoid
    void generate_multi_sphere(int resolution = 5) {
        sub_spheres.clear();
        if (shape == SeedShape::SPHERE) {
            sub_spheres.push_back({{0,0,0}, radius});
            return;
        }
        // Place spheres along the major axis (a-axis) with radii following ellipsoid profile
        double axes[3] = {a, b, c};
        // Sort to find major axis
        int major = 0;
        if (axes[1] > axes[major]) major = 1;
        if (axes[2] > axes[major]) major = 2;

        double L = axes[major];
        int n = std::max(3, resolution);
        for (int i = 0; i < n; i++) {
            double t = -1.0 + 2.0 * (i + 0.5) / n;  // -1 to 1
            vec3 pos(0, 0, 0);
            pos[major] = t * L;
            // Radius at this position: from ellipsoid cross-section
            double frac = 1.0 - t*t;
            if (frac < 0) frac = 0;
            // Use the smaller two axes for the cross-section radius
            double r1 = (major == 0) ? b : ((major == 1) ? a : a);
            double r2 = (major == 0) ? c : ((major == 1) ? c : b);
            double r = std::sqrt(frac) * std::min(r1, r2);
            if (r > 1e-6) {
                sub_spheres.push_back({pos, r});
            }
        }
        if (sub_spheres.empty()) {
            // Fallback: single bounding sphere
            sub_spheres.push_back({{0,0,0}, std::max({a,b,c})});
        }
    }
};

// Per-particle state
struct Particle {
    uint32_t id = 0;
    int seed_type_idx = 0;  // index into SeedType array
    bool active = true;

    // Translational state
    vec3 pos;
    vec3 vel;
    vec3 force;    // accumulated force this step

    // Rotational state
    quat orientation;
    vec3 omega;    // angular velocity (world frame)
    vec3 torque;   // accumulated torque this step

    // Metrics tracking
    double birth_time = 0.0;
    double exit_time = -1.0;
    int wall_hits = 0;
    double cumulative_wall_impulse = 0.0;
    double max_bounce_energy = 0.0;
    vec3 exit_velocity;
    bool has_exited = false;

    // For bounding sphere (multi-sphere ellipsoid)
    double bounding_radius() const; // defined after SeedType known
};

// Particle system: SoA-friendly container
struct ParticleSystem {
    std::vector<Particle> particles;
    std::vector<SeedType> seed_types;
    std::vector<ContactProperties> seed_seed_props;
    std::vector<ContactProperties> seed_wall_props;
    uint32_t next_id = 0;

    ContactProperties default_seed_seed;
    ContactProperties default_seed_wall;

    uint32_t add_particle(const vec3& pos, const vec3& vel, int type_idx, double time) {
        Particle p;
        p.id = next_id++;
        p.seed_type_idx = type_idx;
        p.pos = pos;
        p.vel = vel;
        p.birth_time = time;
        particles.push_back(p);
        return p.id;
    }

    void remove_inactive() {
        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                [](const Particle& p) { return !p.active; }),
            particles.end());
    }

    const SeedType& type_of(const Particle& p) const {
        return seed_types[p.seed_type_idx];
    }

    double bounding_radius(const Particle& p) const {
        const auto& st = type_of(p);
        if (st.shape == SeedShape::SPHERE) return st.radius;
        return std::max({st.a, st.b, st.c});
    }

    const ContactProperties& get_seed_seed_props() const { return default_seed_seed; }
    const ContactProperties& get_seed_wall_props() const { return default_seed_wall; }
};

} // namespace dem
