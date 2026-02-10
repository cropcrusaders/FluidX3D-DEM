#pragma once
// Main DEM simulator: orchestrates particles, collisions, contacts, integration

#include "dem_math.hpp"
#include "particle.hpp"
#include "geometry.hpp"
#include "contact.hpp"
#include "collision.hpp"
#include "airflow.hpp"
#include "injector.hpp"
#include "metrics.hpp"
#include <string>
#include <functional>

namespace dem {

struct SimConfig {
    // Time
    double dt = 1e-5;           // timestep (s)
    double max_time = 2.0;      // max simulation time (s)
    int output_interval = 1000; // steps between console output

    // Physics
    vec3 gravity = {0, 0, -9.81};

    // Aero
    bool enable_airflow = false;
    double drag_Cd = 0.9;          // constant drag coefficient
    bool use_schiller_naumann = false;  // use Re-dependent Cd
    double air_density = 1.225;    // kg/m^3
    std::string airflow_path;      // path to velocity field file
    std::string airflow_format = "binary";  // "binary" or "npy"

    // Stability
    bool auto_dt = true;  // auto-compute dt from stiffness
    double dt_safety = 0.1;  // fraction of critical timestep

    // Determinism
    uint64_t rng_seed = 12345;

    // Output
    std::string output_dir = "output";
    bool save_trajectories = true;
    int max_trajectory_particles = 100;
    double trajectory_sample_interval = 0.001;

    // Tube axis (for lateral velocity computation)
    vec3 tube_axis = {0, 0, -1};
};

class Simulator {
public:
    SimConfig config;
    ParticleSystem particles;
    Geometry geometry;
    ContactModel contact_model;
    CollisionDetector collision_detector;
    AirflowField airflow;
    Injector injector;
    ExitPlane exit_plane;
    MetricsCollector metrics;
    TrajectoryRecorder trajectory_recorder;
    RNG rng;

    // Initialize from config
    void init();

    // Run full simulation
    void run();

    // Single step
    void step();

    // Current state
    double current_time() const { return time_; }
    uint64_t current_step() const { return step_count_; }
    bool is_finished() const;

    // Callbacks
    std::function<void(uint64_t step, double time, int n_particles)> on_step;
    std::function<void(const Particle&)> on_exit;

private:
    double time_ = 0.0;
    uint64_t step_count_ = 0;

    void compute_forces();
    void integrate();
    void check_exits();
    double compute_critical_dt() const;
    void print_status() const;
};

} // namespace dem
