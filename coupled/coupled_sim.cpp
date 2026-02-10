// Two-way coupled LBM-DEM simulation: cotton harvester door-block scenario
//
// FluidX3D (LBM) computes airflow through a shoot with a baffle obstruction.
// DEM computes cotton particle dynamics with drag from the live velocity field.
// Each timestep:
//   1. Read LBM velocity field from GPU
//   2. Update DEM airflow field with live velocities
//   3. Step DEM (computes drag forces on particles)
//   4. Distribute particle reaction forces (-Fd) back to LBM force field
//   5. Step LBM (fluid feels the particle drag reaction)
//
// Build: cd coupled && mkdir build && cd build && cmake .. && make -j$(nproc)
// Requires: OpenCL runtime (GPU driver or CPU fallback)
//
// To use defines.hpp extensions, edit src/defines.hpp BEFORE building:
//   - Uncomment VOLUME_FORCE and FORCE_FIELD (required for two-way coupling)
//   - Comment out BENCHMARK

#include "defines.hpp"
#include "lbm.hpp"
#include "shapes.hpp"
#include "units.hpp"

#include "simulator.hpp"
#include "airflow.hpp"
#include "particle.hpp"
#include "geometry.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// ─── DEM setup for cotton in a door-block scenario ───────────────────────

static void setup_dem(dem::Simulator& sim, double tube_radius_m, double tube_length_m) {
    // Simulation config
    sim.config.dt = 2e-5;
    sim.config.auto_dt = true;
    sim.config.dt_safety = 0.1;
    sim.config.max_time = 5.0;  // 5 seconds to observe blockage
    sim.config.gravity = {0, 0, -9.81};  // gravity in -z
    sim.config.enable_airflow = true;     // will be populated from live LBM data
    sim.config.use_schiller_naumann = true;
    sim.config.air_density = 1.225;
    sim.config.drag_Cd = 0.9;
    sim.config.rng_seed = 42;
    sim.config.output_dir = "output/door_block";
    sim.config.output_interval = 5000;
    sim.config.save_trajectories = true;
    sim.config.max_trajectory_particles = 50;
    sim.config.trajectory_sample_interval = 0.002;
    sim.config.tube_axis = {0, 0, -1};

    // Seed types: cotton boll fragments
    dem::SeedType cotton;
    cotton.name = "cotton_boll";
    cotton.shape = dem::SeedShape::SPHERE;
    cotton.radius = 0.005;  // 10mm diameter
    cotton.mass = 0.002;    // 2g
    sim.particles.seed_types.push_back(cotton);

    // Contact properties: cotton is soft, low bounce
    sim.particles.default_seed_seed = {5000, 1250, 0.25, 0.5, 0.02};
    sim.particles.default_seed_wall = {8000, 2000, 0.15, 0.4, 0.02};

    // Geometry: vertical cylinder with a baffle plate partway down
    // The baffle creates a constriction where cotton can accumulate and block
    dem::AnalyticPrimitive tube;
    tube.type = dem::PrimitiveType::CYLINDER;
    tube.center = {0, 0, tube_length_m * 0.5};  // centered vertically
    tube.axis = {0, 0, -1};
    tube.length = tube_length_m;
    tube.radius = tube_radius_m;
    sim.geometry.add_primitive(tube);

    // Baffle plate: a horizontal plane that blocks part of the tube
    // This represents a "door" or transition that cotton can pile up against
    dem::AnalyticPrimitive baffle;
    baffle.type = dem::PrimitiveType::PLANE;
    baffle.center = {0, 0, tube_length_m * 0.35};  // 35% from bottom
    baffle.plane_normal = {0, 0, 1};
    baffle.plane_offset = tube_length_m * 0.35;
    sim.geometry.add_primitive(baffle);

    // Injector: cotton enters at top of shoot
    sim.injector.config.position = {0, 0, tube_length_m * 0.95};
    sim.injector.config.direction = {0, 0, -1};
    sim.injector.config.seeds_per_second = 30;
    sim.injector.config.max_seeds = 300;
    sim.injector.config.aperture_shape = "circular";
    sim.injector.config.aperture_radius = tube_radius_m * 0.6;
    sim.injector.config.initial_speed = 0.5;
    sim.injector.config.speed_stddev = 0.1;
    sim.injector.config.lateral_speed_stddev = 0.05;
    sim.injector.config.seed_type_idx = 0;

    // Exit plane at bottom
    sim.exit_plane.point = {0, 0, 0.005};
    sim.exit_plane.normal = {0, 0, 1};
}

// ─── Main: two-way coupled LBM-DEM loop ──────────────────────────────────

#ifndef COUPLED_SIM_NO_MAIN  // allow excluding main() for testing

void main_setup(); // forward declaration (required by FluidX3D headers)
void main_setup() {} // stub -- we use our own main loop below

int main() {
    std::cout << "═══════════════════════════════════════════════════\n"
              << "  Two-Way Coupled LBM-DEM: Cotton Door-Block\n"
              << "═══════════════════════════════════════════════════\n\n";

    // ─── Physical parameters ─────────────────────────────────────────
    const float si_D = 0.050f;       // tube diameter (m)
    const float si_L = 0.400f;       // tube length (m)
    const float si_u_avg = 8.0f;     // mean air velocity (m/s)
    const float si_rho = 1.225f;     // air density
    const float si_nu = 1.516e-5f;   // air kinematic viscosity
    const float si_Re = si_u_avg * si_D / si_nu;

    // ─── LBM grid parameters ────────────────────────────────────────
    const float D = 48.0f;           // tube diameter in lattice units
    const float L = D * si_L / si_D; // tube length in lattice units
    const float u_lbm = 0.05f;       // characteristic LBM velocity

    units.set_m_kg_s(D, u_lbm, 1.0f, si_D, si_u_avg, si_rho);
    const float nu = units.nu_from_Re(si_Re, D, u_lbm);
    const float R = 0.5f * D;
    const float f_drive = units.f_from_u_Poiseuille_3D(u_lbm, 1.0f, nu, R);
    const float si_spacing = units.si_x(1.0f);

    const uint Nx_box = to_uint(D) + 4u;
    const uint Ny_box = to_uint(L);
    const uint Nz_box = Nx_box;

    std::cout << "LBM grid: " << Nx_box << "x" << Ny_box << "x" << Nz_box
              << "  Re=" << (int)si_Re << "  spacing=" << si_spacing*1000 << " mm\n";

    // ─── Create LBM ──────────────────────────────────────────────────
    LBM lbm(Nx_box, Ny_box, Nz_box, nu, 0.0f, -f_drive, 0.0f);

    // Geometry: cylindrical tube with baffle plate
    {
        const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
        const float baffle_y = 0.35f * (float)Ny;  // baffle at 35% from bottom
        const float baffle_gap = R * 0.5f;          // half the tube is open

        parallel_for(lbm.get_N(), [&](ulong n) {
            uint x = 0u, y = 0u, z = 0u;
            lbm.coordinates(n, x, y, z);

            // Cylindrical tube wall
            if (!cylinder(x, y, z, lbm.center(), float3(0u, Ny, 0u), R)) {
                lbm.flags[n] = TYPE_S;
            }

            // Baffle plate: solid plane at baffle_y, but with a gap on one side
            // Gap is where x > center (one half of the tube is open)
            if (y >= (uint)baffle_y && y <= (uint)baffle_y + 1u) {
                float dx = (float)x + 0.5f - 0.5f * (float)Nx;
                float dz = (float)z + 0.5f - 0.5f * (float)Nz;
                float r_pos = sqrtf(dx * dx + dz * dz);
                // Only block the half where dx < 0 (leave a gap on dx > 0 side)
                if (r_pos < R && dx < baffle_gap) {
                    lbm.flags[n] = TYPE_S;
                }
            }
        });
    }

    // ─── Create DEM simulator ────────────────────────────────────────
    dem::Simulator dem_sim;
    setup_dem(dem_sim, (double)(si_D * 0.5), (double)si_L);

    // Initialize DEM airflow field to match LBM grid
    const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
    dem::vec3 grid_origin = {
        -0.5 * Nx * (double)si_spacing,
        -0.5 * Ny * (double)si_spacing,
        -0.5 * Nz * (double)si_spacing
    };
    dem_sim.airflow.init_grid(Nx, Ny, Nz, grid_origin, (double)si_spacing);
    dem_sim.config.enable_airflow = true;

    dem_sim.init();

    // ─── Converge LBM flow before coupling ───────────────────────────
    std::cout << "\n[LBM] Converging initial flow field...\n";
    const ulong warmup_steps = std::min(units.t(0.5f), (ulong)50000);
    lbm.run(warmup_steps);
    std::cout << "[LBM] Warmup done (" << warmup_steps << " steps)\n\n";

    // ─── Two-way coupled time loop ───────────────────────────────────
    const float si_dt_dem = (float)dem_sim.config.dt;
    const float si_dt_lbm = units.si_x(1.0f) / units.si_u(1.0f); // ~one LBM timestep in seconds
    // How many LBM steps per DEM step
    const int lbm_steps_per_dem = std::max(1, (int)(si_dt_dem / si_dt_lbm));
    const float u_conv = units.si_u(1.0f);

    std::cout << "Coupling: " << lbm_steps_per_dem << " LBM steps per DEM step\n";
    std::cout << "DEM dt=" << si_dt_dem << " s, LBM dt~" << si_dt_lbm << " s\n\n";

    // Allocate force feedback buffers
    const size_t N_cells = (size_t)Nx * Ny * Nz;
    std::vector<float> fb_fx(N_cells, 0.0f);
    std::vector<float> fb_fy(N_cells, 0.0f);
    std::vector<float> fb_fz(N_cells, 0.0f);

    int print_interval = 500;
    uint64_t coupled_step = 0;

    while (!dem_sim.is_finished()) {
        // 1. Read LBM velocity from GPU to CPU
        lbm.u.read_from_device();

        // 2. Copy LBM velocity (SI units) into DEM airflow field
        //    LBM velocity is in lattice units, convert to SI
        {
            size_t total = N_cells;
            std::vector<float> si_ux(total), si_uy(total), si_uz(total);
            for (size_t n = 0; n < total; n++) {
                si_ux[n] = lbm.u.x[n] * u_conv;
                si_uy[n] = lbm.u.y[n] * u_conv;
                si_uz[n] = lbm.u.z[n] * u_conv;
            }
            dem_sim.airflow.update_from_arrays(si_ux.data(), si_uy.data(), si_uz.data());
        }

        // 3. Step DEM (computes drag forces on particles)
        dem_sim.step_twoway();

        // 4. Map particle drag reaction forces back to LBM force field
#ifdef FORCE_FIELD
        if (!dem_sim.last_drag_forces.empty()) {
            dem_sim.airflow.distribute_forces_to_grid(
                dem_sim.last_drag_forces,
                fb_fx.data(), fb_fy.data(), fb_fz.data()
            );

            // Convert SI force density (N/m^3) to LBM force density units
            const float f_conv = units.f(1.0f);  // SI -> LBM conversion for force/volume
            for (size_t n = 0; n < N_cells; n++) {
                lbm.F.x[n] = fb_fx[n] * f_conv;
                lbm.F.y[n] = fb_fy[n] * f_conv;
                lbm.F.z[n] = fb_fz[n] * f_conv;
            }
            lbm.F.write_to_device();
        }
#endif // FORCE_FIELD

        // 5. Advance LBM
        lbm.run(lbm_steps_per_dem);

        coupled_step++;
        if (coupled_step % print_interval == 0) {
            int active = 0;
            for (auto& p : dem_sim.particles.particles) if (p.active) active++;

            std::cout << "[COUPLED] step=" << coupled_step
                      << " t=" << std::fixed << std::setprecision(4) << dem_sim.current_time()
                      << " active=" << active
                      << " injected=" << dem_sim.injector.total_injected()
                      << " exited=" << dem_sim.metrics.per_seed().size()
                      << " drags=" << dem_sim.last_drag_forces.size()
                      << "\n";
        }
    }

    // ─── Final output ────────────────────────────────────────────────
    std::cout << "\n═══════════════════════════════════════════════════\n"
              << "  Simulation Complete\n"
              << "  Total coupled steps: " << coupled_step << "\n"
              << "  Seeds exited: " << dem_sim.metrics.per_seed().size() << "\n"
              << "═══════════════════════════════════════════════════\n";

    dem_sim.metrics.compute_summary();
    dem_sim.metrics.write_csv(dem_sim.config.output_dir + "/seeds.csv");
    dem_sim.metrics.write_summary(dem_sim.config.output_dir + "/summary.csv");

    // Export final velocity field for analysis
    lbm.write_velocity_binary("output/door_block/final_velocity.bin", si_spacing);

    return 0;
}

#endif // COUPLED_SIM_NO_MAIN
