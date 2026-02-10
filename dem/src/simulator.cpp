// Main DEM simulator: integration loop, force computation, coordination

#include "simulator.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dem {

void Simulator::init() {
    rng = RNG(config.rng_seed);

    // Generate multi-sphere representations
    for (auto& st : particles.seed_types) {
        if (st.shape == SeedShape::MULTI_SPHERE_ELLIPSOID) {
            st.generate_multi_sphere();
        } else {
            st.sub_spheres.clear();
            st.sub_spheres.push_back({{0,0,0}, st.radius});
        }
    }

    // Auto-compute timestep if requested
    if (config.auto_dt) {
        double critical_dt = compute_critical_dt();
        if (critical_dt > 0) {
            config.dt = config.dt_safety * critical_dt;
            std::cout << "[DEM] Auto timestep: dt = " << config.dt
                      << " s (critical = " << critical_dt << " s)\n";
        }
    }

    // Stability warning
    {
        double critical_dt = compute_critical_dt();
        if (critical_dt > 0 && config.dt > critical_dt) {
            std::cerr << "[DEM] WARNING: dt = " << config.dt
                      << " exceeds critical timestep " << critical_dt
                      << " - simulation may be unstable!\n";
        }
    }

    // Load airflow field
    if (config.enable_airflow && !config.airflow_path.empty()) {
        bool ok = false;
        if (config.airflow_format == "npy") {
            ok = airflow.load_npy(config.airflow_path);
        } else {
            ok = airflow.load_binary(config.airflow_path);
        }
        if (ok) {
            airflow.rho_fluid = config.air_density;
            std::cout << "[DEM] Loaded airflow field: "
                      << airflow.nx << "x" << airflow.ny << "x" << airflow.nz
                      << " grid\n";
        } else {
            std::cerr << "[DEM] WARNING: Failed to load airflow field: "
                      << config.airflow_path << "\n";
            config.enable_airflow = false;
        }
    }

    // Setup collision detector
    collision_detector.set_geometry(&geometry);

    // Setup metrics
    metrics.tube_axis = config.tube_axis.normalized();

    // Setup trajectory recorder
    if (config.save_trajectories) {
        trajectory_recorder.init(config.max_trajectory_particles);
        trajectory_recorder.sample_interval = config.trajectory_sample_interval;
    }

    // Init injector
    injector.init(injector.config, rng);

    // Create output directory (recursive)
    {
        std::string path = config.output_dir;
        for (size_t i = 1; i < path.size(); i++) {
            if (path[i] == '/') {
                path[i] = '\0';
                mkdir(path.c_str(), 0755);
                path[i] = '/';
            }
        }
        mkdir(path.c_str(), 0755);
    }

    time_ = 0.0;
    step_count_ = 0;
}

double Simulator::compute_critical_dt() const {
    // dt_crit = sqrt(m_min / k_max) * safety
    double m_min = 1e30;
    for (auto& st : particles.seed_types) {
        m_min = std::min(m_min, st.mass);
    }

    double k_max = std::max(particles.default_seed_seed.kn, particles.default_seed_wall.kn);
    if (k_max <= 0 || m_min <= 0 || m_min > 1e20) return 1e-4;

    return std::sqrt(m_min / k_max);
}

bool Simulator::is_finished() const {
    if (time_ >= config.max_time) return true;
    if (injector.done() && particles.particles.empty()) return true;
    return false;
}

void Simulator::run() {
    std::cout << "[DEM] Starting simulation: max_time=" << config.max_time
              << " dt=" << config.dt << "\n";

    while (!is_finished()) {
        step();

        if (config.output_interval > 0 && step_count_ % config.output_interval == 0) {
            print_status();
        }
    }

    // Final cleanup
    print_status();
    std::cout << "[DEM] Simulation complete. "
              << metrics.per_seed().size() << " seeds exited.\n";

    // Compute and write metrics
    metrics.compute_summary();
    metrics.write_csv(config.output_dir + "/seeds.csv");
    metrics.write_summary(config.output_dir + "/summary.csv");

    if (config.save_trajectories) {
        metrics.write_vtk_trajectories(config.output_dir + "/trajectories.vtk",
                                        trajectory_recorder.positions(),
                                        trajectory_recorder.times());
        metrics.write_trajectory_csv(config.output_dir + "/trajectories.csv",
                                      trajectory_recorder.positions(),
                                      trajectory_recorder.times(),
                                      trajectory_recorder.tracked_ids());
    }
}

void Simulator::step() {
    // 1. Inject new seeds
    injector.inject(particles, time_, config.dt, rng);

    // 2. Detect collisions
    collision_detector.detect(particles);

    // 3. Compute forces
    compute_forces();

    // 4. Integrate
    integrate();

    // 5. Check exits
    check_exits();

    // 6. Record trajectories
    if (config.save_trajectories) {
        trajectory_recorder.sample(particles, time_);
    }

    // 7. Cleanup stale contacts
    contact_model.cleanup_stale_contacts();

    // 8. Remove inactive particles
    particles.remove_inactive();

    time_ += config.dt;
    step_count_++;
}

void Simulator::compute_forces() {
    auto& parts = particles.particles;

    // Reset forces
    for (auto& p : parts) {
        if (!p.active) continue;
        const auto& st = particles.type_of(p);
        p.force = config.gravity * st.mass;
        p.torque = {0,0,0};
    }

    // Seed-seed contact forces
    const auto& ss_props = particles.get_seed_seed_props();
    for (auto& sc : collision_detector.sphere_contacts()) {
        auto& a = parts[sc.idx_a];
        auto& b = parts[sc.idx_b];
        if (!a.active || !b.active) continue;

        const auto& sta = particles.type_of(a);
        const auto& stb = particles.type_of(b);
        double ra = particles.bounding_radius(a);
        double m_eff = (sta.mass * stb.mass) / (sta.mass + stb.mass);

        auto cf = contact_model.compute_sphere_wall(a, ra, sc.contact_point, sc.normal, sc.overlap,
                                                     ss_props, config.dt, m_eff);

        a.force += cf.force;
        b.force -= cf.force;
        a.torque += cf.torque_a;
        b.torque += cf.torque_b;
    }

    // Seed-wall contact forces
    const auto& sw_props = particles.get_seed_wall_props();
    for (auto& wc : collision_detector.wall_contacts()) {
        auto& p = parts[wc.idx];
        if (!p.active) continue;

        const auto& st = particles.type_of(p);
        double r = particles.bounding_radius(p);

        auto cf = contact_model.compute_sphere_wall(p, r, wc.contact_point, wc.normal,
                                                     wc.overlap, sw_props, config.dt, st.mass);

        p.force += cf.force;
        p.torque += cf.torque_a;

        // Metrics: wall hit detection
        double impact_energy = 0.5 * st.mass * p.vel.length2();
        if (cf.normal_impulse > 1e-8) {
            p.wall_hits++;
            p.cumulative_wall_impulse += cf.normal_impulse;
            p.max_bounce_energy = std::max(p.max_bounce_energy, impact_energy);
        }
    }

    // Aero forces (one-way coupling from LBM field)
    if (config.enable_airflow && airflow.is_loaded()) {
        for (auto& p : parts) {
            if (!p.active) continue;
            const auto& st = particles.type_of(p);
            double r = (st.shape == SeedShape::SPHERE) ? st.radius : std::max({st.a, st.b, st.c});

            vec3 fluid_vel = airflow.interpolate_velocity(p.pos);
            vec3 v_rel = fluid_vel - p.vel;
            double v_rel_mag = v_rel.length();

            double Cd = config.drag_Cd;
            if (config.use_schiller_naumann && v_rel_mag > 1e-10) {
                double mu_air = 1.81e-5;  // dynamic viscosity of air
                double Re = config.air_density * v_rel_mag * (2.0 * r) / mu_air;
                Cd = AirflowField::schiller_naumann_Cd(Re);
            }

            vec3 Fd = AirflowField::compute_drag(fluid_vel, p.vel,
                                                   config.air_density, Cd, r);
            p.force += Fd;
        }
    }
}

void Simulator::integrate() {
    double dt = config.dt;

    for (auto& p : particles.particles) {
        if (!p.active) continue;
        const auto& st = particles.type_of(p);

        // Velocity Verlet: v(t+dt/2) = v(t) + a(t)*dt/2
        vec3 acc = p.force / st.mass;
        p.vel += acc * (dt * 0.5);

        // Position update: x(t+dt) = x(t) + v(t+dt/2)*dt
        p.pos += p.vel * dt;

        // Second half of velocity update (forces recomputed next step)
        // For simple leapfrog: v(t+dt) = v(t+dt/2) + a(t)*dt/2
        p.vel += acc * (dt * 0.5);

        // Angular integration
        vec3 inertia = st.inertia();
        vec3 alpha;
        if (inertia.x > 1e-30) alpha.x = p.torque.x / inertia.x;
        if (inertia.y > 1e-30) alpha.y = p.torque.y / inertia.y;
        if (inertia.z > 1e-30) alpha.z = p.torque.z / inertia.z;

        // For non-spherical particles, transform torque to body frame
        if (st.shape != SeedShape::SPHERE) {
            mat3 R = p.orientation.to_matrix();
            mat3 Rt = R.transposed();
            vec3 torque_body = Rt * p.torque;
            if (inertia.x > 1e-30) alpha.x = torque_body.x / inertia.x;
            if (inertia.y > 1e-30) alpha.y = torque_body.y / inertia.y;
            if (inertia.z > 1e-30) alpha.z = torque_body.z / inertia.z;
            // Transform back to world frame
            alpha = R * alpha;
        }

        p.omega += alpha * dt;

        // Quaternion integration
        p.orientation = quat::integrate(p.orientation, p.omega, dt);
    }
}

void Simulator::check_exits() {
    for (auto& p : particles.particles) {
        if (!p.active || p.has_exited) continue;

        if (exit_plane.has_crossed(p.pos)) {
            p.has_exited = true;
            p.exit_time = time_;
            p.exit_velocity = p.vel;

            metrics.record_exit(p, time_);
            if (on_exit) on_exit(p);

            p.active = false;
        }
    }
}

void Simulator::print_status() const {
    int active = 0;
    for (auto& p : particles.particles) if (p.active) active++;

    std::cout << "[DEM] step=" << step_count_
              << " t=" << std::fixed << std::setprecision(4) << time_
              << " active=" << active
              << " injected=" << injector.total_injected()
              << " exited=" << metrics.per_seed().size()
              << " contacts=" << contact_model.active_contacts()
              << "\n";
}

} // namespace dem
