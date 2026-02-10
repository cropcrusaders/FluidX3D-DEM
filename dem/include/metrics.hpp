#pragma once
// Metrics collection and output for DEM simulation

#include "dem_math.hpp"
#include "particle.hpp"
#include <vector>
#include <string>
#include <fstream>

namespace dem {

struct PerSeedMetrics {
    uint32_t id;
    double time_to_exit;
    int wall_hits;
    double cumulative_impulse;
    double max_bounce_energy;
    vec3 exit_velocity;
    double exit_speed;
    double exit_lateral_velocity;
};

struct SimulationSummary {
    int total_seeds;
    int exited_seeds;

    double mean_time_to_exit, sd_time_to_exit;
    double mean_wall_hits, sd_wall_hits;
    double mean_exit_speed, sd_exit_speed;
    double mean_exit_lateral_vel, sd_exit_lateral_vel;
    double mean_cumulative_impulse;
    double mean_max_bounce_energy;

    // Spacing risk proxy
    double spacing_risk_proxy;
};

class MetricsCollector {
public:
    vec3 tube_axis = {0, 0, -1};  // for computing lateral velocity

    void record_exit(const Particle& p, double current_time);
    void compute_summary();

    const std::vector<PerSeedMetrics>& per_seed() const { return per_seed_; }
    const SimulationSummary& summary() const { return summary_; }

    // Output
    bool write_csv(const std::string& path) const;
    bool write_summary(const std::string& path) const;

    // VTK trajectory output
    bool write_vtk_trajectories(const std::string& path,
                                 const std::vector<std::vector<vec3>>& trajectories,
                                 const std::vector<std::vector<double>>& times) const;

    // CSV trajectory output (simpler alternative)
    bool write_trajectory_csv(const std::string& path,
                               const std::vector<std::vector<vec3>>& trajectories,
                               const std::vector<std::vector<double>>& times,
                               const std::vector<uint32_t>& ids) const;

private:
    std::vector<PerSeedMetrics> per_seed_;
    SimulationSummary summary_ = {};
};

// Trajectory recorder: samples particle positions at intervals
class TrajectoryRecorder {
public:
    double sample_interval = 0.001;  // seconds between samples
    int max_particles_to_track = 100;

    void init(int max_track);
    void sample(const ParticleSystem& psys, double current_time);

    const std::vector<std::vector<vec3>>& positions() const { return positions_; }
    const std::vector<std::vector<double>>& times() const { return times_; }
    const std::vector<uint32_t>& tracked_ids() const { return tracked_ids_; }

private:
    std::vector<uint32_t> tracked_ids_;
    std::vector<std::vector<vec3>> positions_;
    std::vector<std::vector<double>> times_;
    double last_sample_time_ = -1e30;
};

} // namespace dem
