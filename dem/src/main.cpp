// DEM Seed Drop Simulator - CLI Entry Point
// Usage: dem_sim <config.json> [--output-dir <dir>] [--quiet]

#include "simulator.hpp"
#include "config_parser.hpp"
#include <iostream>
#include <string>
#include <chrono>

static void print_usage(const char* argv0) {
    std::cout << "DEM Seed Drop-Tube Simulator\n"
              << "Usage: " << argv0 << " <config.json> [options]\n\n"
              << "Options:\n"
              << "  --output-dir <dir>  Override output directory\n"
              << "  --dt <seconds>      Override timestep\n"
              << "  --max-time <sec>    Override max simulation time\n"
              << "  --seed <int>        Override RNG seed\n"
              << "  --quiet             Suppress progress output\n"
              << "  --help              Show this help\n\n"
              << "Output files (in output directory):\n"
              << "  seeds.csv           Per-seed metrics\n"
              << "  summary.csv         Aggregate statistics\n"
              << "  trajectories.vtk    Trajectory data (VTK polydata)\n"
              << "  trajectories.csv    Trajectory data (CSV)\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string config_path;
    std::string output_dir_override;
    double dt_override = -1;
    double max_time_override = -1;
    int64_t seed_override = -1;
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir_override = argv[++i];
        } else if (arg == "--dt" && i + 1 < argc) {
            dt_override = std::stod(argv[++i]);
        } else if (arg == "--max-time" && i + 1 < argc) {
            max_time_override = std::stod(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            seed_override = std::stoll(argv[++i]);
        } else if (arg == "--quiet" || arg == "-q") {
            quiet = true;
        } else if (arg[0] != '-') {
            config_path = arg;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        }
    }

    if (config_path.empty()) {
        std::cerr << "Error: No config file specified\n";
        print_usage(argv[0]);
        return 1;
    }

    // Load config and build simulator
    dem::Simulator sim;

    if (!dem::load_config(config_path, sim)) {
        std::cerr << "Error: Failed to load config from " << config_path << "\n";
        return 1;
    }

    // Apply overrides
    if (!output_dir_override.empty()) sim.config.output_dir = output_dir_override;
    if (dt_override > 0) { sim.config.dt = dt_override; sim.config.auto_dt = false; }
    if (max_time_override > 0) sim.config.max_time = max_time_override;
    if (seed_override >= 0) sim.config.rng_seed = (uint64_t)seed_override;
    if (quiet) sim.config.output_interval = 0;

    // Initialize
    sim.init();

    if (!quiet) {
        std::cout << "═══════════════════════════════════════════════════\n"
                  << "  DEM Seed Drop-Tube Simulator\n"
                  << "═══════════════════════════════════════════════════\n"
                  << "  Config: " << config_path << "\n"
                  << "  Output: " << sim.config.output_dir << "\n"
                  << "  dt = " << sim.config.dt << " s\n"
                  << "  max_time = " << sim.config.max_time << " s\n"
                  << "  max_seeds = " << sim.injector.config.max_seeds << "\n"
                  << "  seed_types = " << sim.particles.seed_types.size() << "\n"
                  << "  geometry: " << sim.geometry.meshes.size() << " meshes, "
                  << sim.geometry.primitives.size() << " primitives\n"
                  << "  airflow: " << (sim.config.enable_airflow ? "enabled" : "disabled") << "\n"
                  << "  RNG seed: " << sim.config.rng_seed << "\n"
                  << "═══════════════════════════════════════════════════\n";
    }

    // Run simulation
    auto t_start = std::chrono::high_resolution_clock::now();
    sim.run();
    auto t_end = std::chrono::high_resolution_clock::now();

    double wall_time = std::chrono::duration<double>(t_end - t_start).count();

    // Print summary
    auto& summary = sim.metrics.summary();
    std::cout << "\n═══════════════════════════════════════════════════\n"
              << "  Results Summary\n"
              << "═══════════════════════════════════════════════════\n"
              << "  Wall clock time: " << wall_time << " s\n"
              << "  Seeds simulated: " << summary.total_seeds << "\n"
              << "  Seeds exited:    " << summary.exited_seeds << "\n"
              << "  ─────────────────────────────────────────────────\n"
              << "  Time to exit:    " << summary.mean_time_to_exit
              << " +/- " << summary.sd_time_to_exit << " s\n"
              << "  Wall hits:       " << summary.mean_wall_hits
              << " +/- " << summary.sd_wall_hits << "\n"
              << "  Exit speed:      " << summary.mean_exit_speed
              << " +/- " << summary.sd_exit_speed << " m/s\n"
              << "  Exit lateral v:  " << summary.mean_exit_lateral_vel
              << " +/- " << summary.sd_exit_lateral_vel << " m/s\n"
              << "  Avg impulse:     " << summary.mean_cumulative_impulse << " N*s\n"
              << "  Avg max bounce:  " << summary.mean_max_bounce_energy << " J\n"
              << "  ─────────────────────────────────────────────────\n"
              << "  SPACING RISK PROXY: " << summary.spacing_risk_proxy << "\n"
              << "  (lower = better singulation)\n"
              << "═══════════════════════════════════════════════════\n"
              << "\n  Output files:\n"
              << "    " << sim.config.output_dir << "/seeds.csv\n"
              << "    " << sim.config.output_dir << "/summary.csv\n";

    if (sim.config.save_trajectories) {
        std::cout << "    " << sim.config.output_dir << "/trajectories.vtk\n"
                  << "    " << sim.config.output_dir << "/trajectories.csv\n";
    }

    std::cout << "\n";
    return 0;
}
