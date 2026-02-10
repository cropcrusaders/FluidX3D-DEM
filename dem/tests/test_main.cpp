// DEM Validation Tests
// Tests: free-fall, restitution, friction, energy conservation, math, config

#include "dem_math.hpp"
#include "particle.hpp"
#include "geometry.hpp"
#include "contact.hpp"
#include "collision.hpp"
#include "airflow.hpp"
#include "injector.hpp"
#include "metrics.hpp"
#include "simulator.hpp"
#include "config_parser.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        std::cerr << "  FAIL: " << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define CHECK_NEAR(a, b, tol, msg) do { \
    double _a = (a), _b = (b), _t = (tol); \
    if (std::abs(_a - _b) > _t) { \
        std::cerr << "  FAIL: " << msg << " (got " << _a << ", expected " << _b \
                  << ", tol " << _t << ")\n"; \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

// ─── Test 1: Vector math ────────────────────────────────────────────────────

void test_vec3_math() {
    std::cout << "Test: vec3 math\n";
    using namespace dem;

    vec3 a(1, 2, 3), b(4, 5, 6);
    vec3 c = a + b;
    CHECK_NEAR(c.x, 5, 1e-10, "vec3 add x");
    CHECK_NEAR(c.y, 7, 1e-10, "vec3 add y");
    CHECK_NEAR(c.z, 9, 1e-10, "vec3 add z");

    CHECK_NEAR(a.dot(b), 32, 1e-10, "vec3 dot");
    vec3 cr = a.cross(b);
    CHECK_NEAR(cr.x, -3, 1e-10, "vec3 cross x");
    CHECK_NEAR(cr.y, 6, 1e-10, "vec3 cross y");
    CHECK_NEAR(cr.z, -3, 1e-10, "vec3 cross z");

    CHECK_NEAR(vec3(3,4,0).length(), 5, 1e-10, "vec3 length");
    vec3 n = vec3(0,0,5).normalized();
    CHECK_NEAR(n.z, 1.0, 1e-10, "vec3 normalized");
}

// ─── Test 2: Quaternion ─────────────────────────────────────────────────────

void test_quaternion() {
    std::cout << "Test: quaternion\n";
    using namespace dem;

    quat q(1, 0, 0, 0);  // identity
    vec3 v(1, 0, 0);
    vec3 r = q.rotate(v);
    CHECK_NEAR(r.x, 1, 1e-10, "quat identity rotate x");
    CHECK_NEAR(r.y, 0, 1e-10, "quat identity rotate y");

    // 90-degree rotation around Z
    double a = M_PI / 2.0;
    quat rz(std::cos(a/2), 0, 0, std::sin(a/2));
    vec3 rotated = rz.rotate(vec3(1, 0, 0));
    CHECK_NEAR(rotated.x, 0, 1e-10, "quat 90z rotate x");
    CHECK_NEAR(rotated.y, 1, 1e-10, "quat 90z rotate y");

    // Normalization
    quat q2(2, 0, 0, 0);
    quat n = q2.normalized();
    CHECK_NEAR(n.w, 1.0, 1e-10, "quat normalize");
}

// ─── Test 3: Free-fall (no contacts) ────────────────────────────────────────

void test_free_fall() {
    std::cout << "Test: free-fall trajectory\n";
    using namespace dem;

    // Single particle falling under gravity, no contacts
    // z(t) = z0 + v0*t - 0.5*g*t^2
    double g = 9.81;
    double z0 = 1.0;
    double v0 = 0.0;
    double t_total = 0.4;  // ~0.78m fall

    Simulator sim;
    sim.config.dt = 1e-5;
    sim.config.auto_dt = false;
    sim.config.max_time = t_total;
    sim.config.gravity = {0, 0, -g};
    sim.config.output_interval = 0;
    sim.config.save_trajectories = false;
    sim.config.output_dir = "/tmp/dem_test_freefall";

    SeedType st;
    st.radius = 0.004;
    st.mass = 0.001;
    sim.particles.seed_types.push_back(st);
    sim.particles.default_seed_seed = {1e4, 2500, 0.5, 0.5, 0.01};
    sim.particles.default_seed_wall = {1e4, 2500, 0.4, 0.4, 0.01};

    // No geometry, no exit plane (set exit very far below)
    sim.exit_plane.point = {0, 0, -100};
    sim.exit_plane.normal = {0, 0, 1};

    // Disable injector, manually add particle
    sim.injector.config.max_seeds = 0;
    sim.particles.add_particle({0, 0, z0}, {0, 0, v0}, 0, 0.0);

    sim.init();

    // Run
    while (!sim.is_finished()) {
        sim.step();
    }

    // Check final position
    double t = sim.current_time();
    double z_expected = z0 + v0 * t - 0.5 * g * t * t;
    double z_actual = sim.particles.particles[0].pos.z;

    CHECK_NEAR(z_actual, z_expected, 0.001, "free-fall z position");

    // Check velocity
    double vz_expected = v0 - g * t;
    double vz_actual = sim.particles.particles[0].vel.z;
    CHECK_NEAR(vz_actual, vz_expected, 0.01, "free-fall z velocity");
}

// ─── Test 4: Bounce on flat plane (restitution) ─────────────────────────────

void test_restitution() {
    std::cout << "Test: bounce restitution\n";
    using namespace dem;

    double e_target = 0.6;
    double drop_height = 0.1;  // m
    double g = 9.81;

    Simulator sim;
    sim.config.dt = 1e-6;  // small dt for accurate bounce
    sim.config.auto_dt = false;
    sim.config.max_time = 0.5;
    sim.config.gravity = {0, 0, -g};
    sim.config.output_interval = 0;
    sim.config.save_trajectories = false;
    sim.config.output_dir = "/tmp/dem_test_restitution";

    SeedType st;
    st.radius = 0.004;
    st.mass = 0.001;
    sim.particles.seed_types.push_back(st);

    // Wall properties with target restitution
    ContactProperties wall_props;
    wall_props.kn = 1e5;
    wall_props.kt = 2500;
    wall_props.restitution = e_target;
    wall_props.friction_s = 0.0;  // no friction for this test
    wall_props.friction_r = 0.0;
    sim.particles.default_seed_wall = wall_props;
    sim.particles.default_seed_seed = wall_props;

    // Floor plane at z=0
    AnalyticPrimitive floor;
    floor.type = PrimitiveType::PLANE;
    floor.center = {0, 0, 0};
    floor.plane_normal = {0, 0, 1};
    floor.plane_offset = 0;
    sim.geometry.add_primitive(floor);

    // No exit, no injector
    sim.exit_plane.point = {0, 0, -100};
    sim.exit_plane.normal = {0, 0, 1};
    sim.injector.config.max_seeds = 0;

    // Drop particle from height
    sim.particles.add_particle({0, 0, drop_height + st.radius}, {0, 0, 0}, 0, 0.0);

    sim.init();

    // Track velocity: find min vz (max approach speed), then max vz after bounce
    double min_vz = 0;   // minimum (most negative) velocity before first contact
    double max_vz_after = 0;  // maximum (most positive) velocity after first zero-crossing
    bool was_negative = false;
    bool crossed_zero = false;

    while (!sim.is_finished() && sim.current_time() < 0.5) {
        sim.step();
        auto& p = sim.particles.particles[0];

        if (!crossed_zero) {
            if (p.vel.z < min_vz) min_vz = p.vel.z;
            if (p.vel.z < -0.1) was_negative = true;
            if (was_negative && p.vel.z > 0) {
                crossed_zero = true;
                max_vz_after = p.vel.z;
            }
        } else {
            // Track max upward velocity after first bounce
            if (p.vel.z > max_vz_after) max_vz_after = p.vel.z;
            // Stop once it starts coming back down
            if (p.vel.z < max_vz_after - 0.01 && max_vz_after > 0.1) break;
        }
    }

    CHECK(crossed_zero, "bounce detected");
    if (crossed_zero) {
        double v_impact = std::abs(min_vz);
        double v_rebound = max_vz_after;
        double e_measured = v_rebound / v_impact;
        CHECK_NEAR(e_measured, e_target, 0.2, "restitution coefficient");
    }
}

// ─── Test 5: Inclined plane friction ────────────────────────────────────────

void test_friction_incline() {
    std::cout << "Test: inclined plane friction\n";
    using namespace dem;

    // A sphere on an inclined plane should slide if tan(angle) > mu
    // and stay put if tan(angle) < mu
    double mu = 0.5;
    double angle_slide = std::atan(mu) + 0.1;  // should slide
    double g = 9.81;

    Simulator sim;
    sim.config.dt = 1e-5;
    sim.config.auto_dt = false;
    sim.config.max_time = 0.1;
    sim.config.gravity = {g * std::sin(angle_slide), 0, -g * std::cos(angle_slide)};
    sim.config.output_interval = 0;
    sim.config.save_trajectories = false;
    sim.config.output_dir = "/tmp/dem_test_friction";

    SeedType st;
    st.radius = 0.004;
    st.mass = 0.001;
    sim.particles.seed_types.push_back(st);

    ContactProperties wall_props;
    wall_props.kn = 1e5;
    wall_props.kt = 5e4;
    wall_props.restitution = 0.0;
    wall_props.friction_s = mu;
    wall_props.friction_r = 0.0;
    sim.particles.default_seed_wall = wall_props;
    sim.particles.default_seed_seed = wall_props;

    AnalyticPrimitive floor;
    floor.type = PrimitiveType::PLANE;
    floor.center = {0, 0, 0};
    floor.plane_normal = {0, 0, 1};
    floor.plane_offset = 0;
    sim.geometry.add_primitive(floor);

    sim.exit_plane.point = {0, 0, -100};
    sim.exit_plane.normal = {0, 0, 1};
    sim.injector.config.max_seeds = 0;

    sim.particles.add_particle({0, 0, st.radius + 0.001}, {0, 0, 0}, 0, 0.0);

    sim.init();

    while (!sim.is_finished()) {
        sim.step();
    }

    auto& p = sim.particles.particles[0];
    // Should have moved in the x direction (sliding downhill)
    CHECK(p.pos.x > 0.001, "particle slides on steep incline");
    CHECK(p.vel.x > 0.0, "particle has positive x velocity on steep incline");
}

// ─── Test 6: Energy conservation check ──────────────────────────────────────

void test_energy_conservation() {
    std::cout << "Test: energy conservation (no dissipation)\n";
    using namespace dem;

    // Two spheres colliding head-on with e=1 (elastic)
    // Total kinetic energy should be conserved
    Simulator sim;
    sim.config.dt = 1e-6;
    sim.config.auto_dt = false;
    sim.config.max_time = 0.01;
    sim.config.gravity = {0, 0, 0};  // no gravity
    sim.config.output_interval = 0;
    sim.config.save_trajectories = false;
    sim.config.output_dir = "/tmp/dem_test_energy";

    SeedType st;
    st.radius = 0.005;
    st.mass = 0.001;
    sim.particles.seed_types.push_back(st);

    ContactProperties props;
    props.kn = 1e5;
    props.kt = 2500;
    props.restitution = 1.0;  // perfectly elastic
    props.friction_s = 0.0;
    props.friction_r = 0.0;
    sim.particles.default_seed_seed = props;
    sim.particles.default_seed_wall = props;

    sim.exit_plane.point = {0, 0, -100};
    sim.exit_plane.normal = {0, 0, 1};
    sim.injector.config.max_seeds = 0;

    // Two particles approaching each other
    double v0 = 1.0;
    sim.particles.add_particle({-0.02, 0, 0}, {v0, 0, 0}, 0, 0.0);
    sim.particles.add_particle({ 0.02, 0, 0}, {-v0, 0, 0}, 0, 0.0);

    sim.init();

    double KE_initial = 0.5 * st.mass * v0 * v0 * 2;

    while (!sim.is_finished()) {
        sim.step();
    }

    double KE_final = 0;
    for (auto& p : sim.particles.particles) {
        KE_final += 0.5 * st.mass * p.vel.length2();
    }

    // Energy should be approximately conserved (within 10% - some numerical loss expected)
    CHECK_NEAR(KE_final, KE_initial, KE_initial * 0.1, "elastic collision energy conservation");
}

// ─── Test 7: RNG determinism ────────────────────────────────────────────────

void test_rng_determinism() {
    std::cout << "Test: RNG determinism\n";
    using namespace dem;

    RNG rng1(42), rng2(42);
    bool all_same = true;
    for (int i = 0; i < 1000; i++) {
        if (rng1.uniform() != rng2.uniform()) {
            all_same = false;
            break;
        }
    }
    CHECK(all_same, "RNG with same seed produces same sequence");

    RNG rng3(42), rng4(99);
    bool any_different = false;
    for (int i = 0; i < 100; i++) {
        if (rng3.uniform() != rng4.uniform()) {
            any_different = true;
            break;
        }
    }
    CHECK(any_different, "RNG with different seeds produces different sequences");
}

// ─── Test 8: JSON parsing ───────────────────────────────────────────────────

void test_json_parser() {
    std::cout << "Test: JSON parser\n";
    using namespace dem;

    std::string json = R"({
        "simulation": {
            "dt": 1e-5,
            "max_time": 2.0,
            "gravity": [0, 0, -9.81],
            "rng_seed": 12345
        },
        "seed_types": [
            {
                "name": "cotton",
                "shape": "sphere",
                "radius": 0.004,
                "mass": 0.001
            }
        ],
        "injector": {
            "position": [0, 0, 0.3],
            "max_seeds": 100
        }
    })";

    Simulator sim;
    bool ok = parse_config(json, sim);
    CHECK(ok, "JSON parse succeeds");
    CHECK_NEAR(sim.config.max_time, 2.0, 1e-10, "JSON max_time parsed");
    CHECK_NEAR(sim.config.gravity.z, -9.81, 1e-10, "JSON gravity parsed");
    CHECK(sim.particles.seed_types.size() == 1, "JSON seed type count");
    CHECK(sim.particles.seed_types[0].name == "cotton", "JSON seed type name");
    CHECK(sim.injector.config.max_seeds == 100, "JSON injector max_seeds");
}

// ─── Test 9: AABB operations ────────────────────────────────────────────────

void test_aabb() {
    std::cout << "Test: AABB operations\n";
    using namespace dem;

    AABB a({0,0,0}, {1,1,1});
    AABB b({0.5,0.5,0.5}, {1.5,1.5,1.5});
    AABB c({2,2,2}, {3,3,3});

    CHECK(a.overlaps(b), "overlapping AABBs");
    CHECK(!a.overlaps(c), "non-overlapping AABBs");
    CHECK(a.contains(vec3(0.5, 0.5, 0.5)), "AABB contains point");
    CHECK(!a.contains(vec3(1.5, 0.5, 0.5)), "AABB does not contain point");
}

// ─── Test 10: Closest point on triangle ─────────────────────────────────────

void test_triangle_closest() {
    std::cout << "Test: closest point on triangle\n";
    using namespace dem;

    Triangle tri;
    tri.v0 = {0, 0, 0};
    tri.v1 = {1, 0, 0};
    tri.v2 = {0, 1, 0};

    // Point above triangle center
    auto r = closest_point_on_triangle({0.25, 0.25, 1.0}, tri);
    CHECK(r.valid, "triangle closest point valid");
    CHECK_NEAR(r.closest_point.z, 0, 1e-10, "closest point on triangle z");
    CHECK_NEAR(r.distance, 1.0, 1e-10, "distance to triangle");

    // Point near vertex
    auto r2 = closest_point_on_triangle({-0.1, -0.1, 0}, tri);
    CHECK(r2.valid, "triangle closest near vertex valid");
    CHECK_NEAR(r2.closest_point.x, 0, 1e-10, "closest to vertex x");
    CHECK_NEAR(r2.closest_point.y, 0, 1e-10, "closest to vertex y");
}

// ─── Test 11: Analytic cylinder collision ───────────────────────────────────

void test_cylinder_collision() {
    std::cout << "Test: analytic cylinder collision\n";
    using namespace dem;

    AnalyticPrimitive cyl;
    cyl.type = PrimitiveType::CYLINDER;
    cyl.center = {0, 0, 0};
    cyl.axis = {0, 0, 1};
    cyl.length = 1.0;
    cyl.radius = 0.05;

    // Particle overlapping wall of cylinder: radial=0.047 + radius=0.004 = 0.051 > 0.05
    auto r = cyl.collide_sphere({0.047, 0, 0.5}, 0.004);
    CHECK(r.hit, "cylinder collision detected");
    CHECK(r.penetration > 0, "positive penetration");

    // Particle well inside (no contact with wall)
    auto r2 = cyl.collide_sphere({0.01, 0, 0.5}, 0.004);
    CHECK(!r2.hit, "no collision when far from wall");

    // Particle outside cylinder length
    auto r3 = cyl.collide_sphere({0.04, 0, 1.5}, 0.004);
    CHECK(!r3.hit, "no collision outside cylinder length");
}

// ─── Test 12: Drag force computation ────────────────────────────────────────

void test_drag_force() {
    std::cout << "Test: drag force computation\n";
    using namespace dem;

    vec3 fluid_vel = {0, 0, -5.0};  // downward airflow
    vec3 particle_vel = {0, 0, -1.0};  // slower particle
    double rho = 1.225;
    double Cd = 0.9;
    double r = 0.004;

    vec3 Fd = AirflowField::compute_drag(fluid_vel, particle_vel, rho, Cd, r);

    // Drag should be in the -z direction (pulling particle downward with air)
    CHECK(Fd.z < 0, "drag force in direction of relative flow");
    CHECK(std::abs(Fd.x) < 1e-15, "no lateral drag");
    CHECK(std::abs(Fd.y) < 1e-15, "no lateral drag");

    // Check magnitude: Fd = 0.5 * rho * Cd * A * |v_rel|^2
    double A = M_PI * r * r;
    double v_rel = 4.0;
    double expected_Fd = 0.5 * rho * Cd * A * v_rel * v_rel;
    CHECK_NEAR(std::abs(Fd.z), expected_Fd, 1e-10, "drag force magnitude");
}

// ─── Test 13: Schiller-Naumann Cd ───────────────────────────────────────────

void test_schiller_naumann() {
    std::cout << "Test: Schiller-Naumann drag coefficient\n";
    using namespace dem;

    // Stokes regime (Re << 1): Cd ~ 24/Re
    double Cd_low = AirflowField::schiller_naumann_Cd(0.1);
    CHECK_NEAR(Cd_low, 24.0/0.1 * (1 + 0.15 * std::pow(0.1, 0.687)), 1.0, "S-N low Re");

    // Newton regime (Re > 1000): Cd = 0.44
    double Cd_high = AirflowField::schiller_naumann_Cd(5000);
    CHECK_NEAR(Cd_high, 0.44, 1e-10, "S-N high Re");
}

// ─── Test 14: Seed type inertia ─────────────────────────────────────────────

void test_inertia() {
    std::cout << "Test: seed inertia\n";
    using namespace dem;

    SeedType sphere;
    sphere.shape = SeedShape::SPHERE;
    sphere.radius = 0.004;
    sphere.mass = 0.001;
    vec3 I = sphere.inertia();
    double expected = 0.4 * 0.001 * 0.004 * 0.004;
    CHECK_NEAR(I.x, expected, 1e-12, "sphere Ixx");
    CHECK_NEAR(I.y, expected, 1e-12, "sphere Iyy");
    CHECK_NEAR(I.z, expected, 1e-12, "sphere Izz");

    SeedType ellipsoid;
    ellipsoid.shape = SeedShape::MULTI_SPHERE_ELLIPSOID;
    ellipsoid.mass = 0.001;
    ellipsoid.a = 0.006; ellipsoid.b = 0.004; ellipsoid.c = 0.003;
    vec3 Ie = ellipsoid.inertia();
    CHECK_NEAR(Ie.x, 0.2 * 0.001 * (0.004*0.004 + 0.003*0.003), 1e-12, "ellipsoid Ixx");
}

// ─── Test: Two-way coupling (init_grid, update_from_arrays, force distribution)

void test_twoway_coupling() {
    std::cout << "Test: two-way coupling (airflow grid update + force feedback)\n";
    using namespace dem;

    // Create a small 4x4x4 airflow grid
    AirflowField field;
    field.init_grid(4, 4, 4, {0, 0, 0}, 0.01);  // 10mm cells, origin at (0,0,0)
    CHECK(field.is_loaded(), "init_grid creates loaded field");
    CHECK(field.nx == 4 && field.ny == 4 && field.nz == 4, "init_grid dimensions");

    // Fill with uniform velocity (1 m/s in z-direction)
    size_t total = 4 * 4 * 4;
    std::vector<float> ux_data(total, 0.0f);
    std::vector<float> uy_data(total, 0.0f);
    std::vector<float> uz_data(total, 2.0f);
    field.update_from_arrays(ux_data.data(), uy_data.data(), uz_data.data());

    // Interpolate at center of grid: should get (0, 0, 2)
    vec3 center = {0.015, 0.015, 0.015};  // center of 4x4x4 grid with 10mm spacing
    vec3 vel = field.interpolate_velocity(center);
    CHECK_NEAR(vel.x, 0.0, 1e-6, "updated field vel.x");
    CHECK_NEAR(vel.z, 2.0, 1e-6, "updated field vel.z");

    // Test force distribution: single particle at grid center
    AirflowField::ParticleDrag pd;
    pd.pos = center;
    pd.drag = {0, 0, 1.0};  // 1N drag in z-direction
    pd.radius = 0.005;

    std::vector<float> fx(total, 0.0f), fy(total, 0.0f), fz(total, 0.0f);
    field.distribute_forces_to_grid({pd}, fx.data(), fy.data(), fz.data());

    // Reaction force is -drag, so fz should be negative (reaction = -1N in z)
    // Sum of all fz entries * cell_volume should equal -1.0 N
    double cell_vol = 0.01 * 0.01 * 0.01;  // 10mm^3
    double total_fz = 0.0;
    for (size_t i = 0; i < total; i++) total_fz += fz[i] * cell_vol;
    CHECK_NEAR(total_fz, -1.0, 1e-6, "force feedback conserves total force");

    // fx and fy should sum to zero (drag was only in z)
    double total_fx = 0.0, total_fy = 0.0;
    for (size_t i = 0; i < total; i++) {
        total_fx += fx[i] * cell_vol;
        total_fy += fy[i] * cell_vol;
    }
    CHECK_NEAR(total_fx, 0.0, 1e-10, "no spurious fx feedback");
    CHECK_NEAR(total_fy, 0.0, 1e-10, "no spurious fy feedback");
}

// ─── Main ───────────────────────────────────────────────────────────────────

int main() {
    std::cout << "═══════════════════════════════════════════════════\n"
              << "  DEM Seed Simulator - Validation Tests\n"
              << "═══════════════════════════════════════════════════\n\n";

    test_vec3_math();
    test_quaternion();
    test_aabb();
    test_triangle_closest();
    test_rng_determinism();
    test_json_parser();
    test_inertia();
    test_drag_force();
    test_schiller_naumann();
    test_cylinder_collision();
    test_free_fall();
    test_restitution();
    test_friction_incline();
    test_energy_conservation();
    test_twoway_coupling();

    std::cout << "\n═══════════════════════════════════════════════════\n"
              << "  Results: " << tests_passed << " passed, "
              << tests_failed << " failed\n"
              << "═══════════════════════════════════════════════════\n";

    return tests_failed > 0 ? 1 : 0;
}
