// Contact model implementation: linear spring-dashpot with Coulomb friction

#include "contact.hpp"
#include <cmath>
#include <algorithm>

namespace dem {

double ContactModel::damping_from_restitution(double e, double m_eff, double kn) {
    // cn = -2 * ln(e) * sqrt(m_eff * kn) / sqrt(pi^2 + ln(e)^2)
    if (e <= 0.0) return 2.0 * std::sqrt(m_eff * kn);  // critically damped
    if (e >= 1.0) return 0.0;  // no damping
    double ln_e = std::log(e);
    return -2.0 * ln_e * std::sqrt(m_eff * kn) / std::sqrt(M_PI * M_PI + ln_e * ln_e);
}

ContactState& ContactModel::get_or_create(uint32_t id_a, uint32_t id_b) {
    // Ensure ordered key
    if (id_a > id_b && id_b != WALL_ID) std::swap(id_a, id_b);
    ContactKey key{id_a, id_b};
    auto& state = contact_history[key];
    state.active = true;
    state.age = 0;
    return state;
}

ContactForce ContactModel::compute_contact(
    double overlap, const vec3& normal,
    const vec3& vel_rel,
    const vec3& omega_a, double ra,
    const vec3& omega_b, double rb,
    double m_eff,
    const ContactProperties& props,
    ContactState& state,
    double dt)
{
    ContactForce result;
    if (overlap <= 0) return result;

    double kn = props.kn;
    double kt = props.kt;
    double cn = damping_from_restitution(props.restitution, m_eff, kn);
    double ct = cn * 0.5;  // tangential damping ~ half normal damping

    // Normal component of relative velocity
    double vn = dot(vel_rel, normal);
    vec3 vel_n = normal * vn;
    vec3 vel_t = vel_rel - vel_n;

    // Add rotational contribution to tangential velocity at contact
    double r_eff = (rb > 0) ? (ra * rb) / (ra + rb) : ra;
    vec3 omega_contact = cross(omega_a, normal * ra);
    if (rb > 0) omega_contact += cross(omega_b, normal * (-rb));
    vel_t = vel_t + omega_contact;

    // Normal force: spring + damping
    // vn = dot(vel_rel, normal): negative when approaching, positive when separating
    // Damping opposes relative motion:
    //   approach (vn<0): -cn*vn > 0 → adds to repulsive force → opposes approach
    //   separation (vn>0): -cn*vn < 0 → reduces repulsive force → opposes separation
    double fn_spring = kn * overlap;
    double fn_total = fn_spring - cn * vn;
    if (fn_total < 0) fn_total = 0;  // no tensile contact

    vec3 Fn = normal * fn_total;

    // Tangential force with history
    // Update tangential displacement
    state.tangential_disp += vel_t * dt;

    // Project tangential displacement onto contact plane
    double proj = dot(state.tangential_disp, normal);
    state.tangential_disp = state.tangential_disp - normal * proj;

    // Tangential spring + damping
    vec3 Ft = state.tangential_disp * (-kt) - vel_t * ct;

    // Coulomb friction limit
    double ft_mag = Ft.length();
    double ft_max = props.friction_s * fn_total;
    if (ft_mag > ft_max && ft_mag > 1e-30) {
        Ft = Ft * (ft_max / ft_mag);
        // Reset tangential displacement to match sliding
        state.tangential_disp = Ft / (-kt);
    }

    result.force = Fn + Ft;

    // Torque from tangential force
    result.torque_a = cross(normal * (-ra), Ft);
    result.torque_b = cross(normal * rb, Ft);

    // Rolling resistance torque
    if (props.friction_r > 0) {
        vec3 omega_rel = omega_a - omega_b;
        double omega_mag = omega_rel.length();
        if (omega_mag > 1e-12) {
            vec3 omega_dir = omega_rel / omega_mag;
            double Tr_mag = props.friction_r * fn_total * r_eff;
            vec3 Tr = omega_dir * (-Tr_mag);
            result.torque_a = result.torque_a + Tr;
            result.torque_b = result.torque_b - Tr;
        }
    }

    result.normal_impulse = fn_total * dt;
    return result;
}

ContactForce ContactModel::compute_sphere_sphere(
    const Particle& a, double ra,
    const Particle& b, double rb,
    const ContactProperties& props,
    double dt)
{
    vec3 diff = b.pos - a.pos;
    double dist = diff.length();
    double overlap = (ra + rb) - dist;

    if (overlap <= 0 || dist < 1e-30) return {};

    vec3 normal = diff / dist;  // from a to b
    vec3 vel_rel = a.vel - b.vel;  // relative velocity of a w.r.t. b

    double m_a = 1.0, m_b = 1.0;  // mass set externally, use 1.0 for now
    double m_eff = (m_a * m_b) / (m_a + m_b);

    auto& state = get_or_create(a.id, b.id);
    return compute_contact(overlap, normal, vel_rel, a.omega, ra, b.omega, rb, m_eff, props, state, dt);
}

ContactForce ContactModel::compute_sphere_wall(
    const Particle& p, double radius,
    const vec3& contact_point,
    const vec3& wall_normal,
    double penetration,
    const ContactProperties& props,
    double dt,
    double particle_mass)
{
    if (penetration <= 0) return {};

    vec3 vel_rel = p.vel;  // wall is static
    double m_eff = particle_mass;  // effective mass = particle mass (wall has infinite mass)

    vec3 zero_omega(0,0,0);
    auto& state = get_or_create(p.id, WALL_ID);
    return compute_contact(penetration, wall_normal, vel_rel, p.omega, radius, zero_omega, 0, m_eff, props, state, dt);
}

void ContactModel::cleanup_stale_contacts(int max_age) {
    for (auto it = contact_history.begin(); it != contact_history.end(); ) {
        it->second.age++;
        if (it->second.age > max_age) {
            it = contact_history.erase(it);
        } else {
            ++it;
        }
    }
}

void ContactModel::reset() {
    contact_history.clear();
}

} // namespace dem
