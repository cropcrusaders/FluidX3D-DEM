#pragma once
// Contact model: linear spring-dashpot with Coulomb friction
// Tracks tangential displacement history per contact pair

#include "dem_math.hpp"
#include "particle.hpp"
#include <unordered_map>
#include <cstdint>

namespace dem {

// Contact pair key (ordered pair of IDs; wall contacts use WALL_ID)
static constexpr uint32_t WALL_ID = 0xFFFFFFFF;

struct ContactKey {
    uint32_t id_a, id_b;
    bool operator==(const ContactKey& o) const { return id_a == o.id_a && id_b == o.id_b; }
};

struct ContactKeyHash {
    size_t operator()(const ContactKey& k) const {
        uint64_t combined = ((uint64_t)k.id_a << 32) | k.id_b;
        // FNV-1a
        combined ^= combined >> 33;
        combined *= 0xff51afd7ed558ccd;
        combined ^= combined >> 33;
        return (size_t)combined;
    }
};

struct ContactState {
    vec3 tangential_disp = {0,0,0};  // accumulated tangential displacement
    bool active = false;
    int age = 0;  // steps since created
};

// Contact force result
struct ContactForce {
    vec3 force;       // force on particle A (reaction on B or wall)
    vec3 torque_a;    // torque on particle A
    vec3 torque_b;    // torque on particle B (zero for wall)
    double normal_impulse = 0.0;  // |Fn| * dt for metrics
};

class ContactModel {
public:
    // Compute contact force between two spheres
    ContactForce compute_sphere_sphere(
        const Particle& a, double ra,
        const Particle& b, double rb,
        const ContactProperties& props,
        double dt);

    // Compute contact force between sphere and wall
    ContactForce compute_sphere_wall(
        const Particle& p, double radius,
        const vec3& contact_point,
        const vec3& wall_normal,
        double penetration,
        const ContactProperties& props,
        double dt,
        double particle_mass = 1.0);

    // Age out stale contacts (call once per step)
    void cleanup_stale_contacts(int max_age = 3);

    // Clear all history
    void reset();

    size_t active_contacts() const { return contact_history.size(); }

private:
    std::unordered_map<ContactKey, ContactState, ContactKeyHash> contact_history;

    ContactState& get_or_create(uint32_t id_a, uint32_t id_b);
    ContactForce compute_contact(
        double overlap, const vec3& normal,
        const vec3& vel_rel,       // relative velocity at contact point
        const vec3& omega_a, double ra,
        const vec3& omega_b, double rb,
        double m_eff,
        const ContactProperties& props,
        ContactState& state,
        double dt);

    // Damping coefficient from restitution
    static double damping_from_restitution(double e, double m_eff, double kn);
};

} // namespace dem
