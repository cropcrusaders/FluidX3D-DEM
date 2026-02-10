// Metrics collection and output

#include "metrics.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>

namespace dem {

void MetricsCollector::record_exit(const Particle& p, double current_time) {
    PerSeedMetrics m;
    m.id = p.id;
    m.time_to_exit = current_time - p.birth_time;
    m.wall_hits = p.wall_hits;
    m.cumulative_impulse = p.cumulative_wall_impulse;
    m.max_bounce_energy = p.max_bounce_energy;
    m.exit_velocity = p.vel;
    m.exit_speed = p.vel.length();

    // Lateral velocity: component perpendicular to tube axis
    double along = dot(p.vel, tube_axis);
    vec3 lateral = p.vel - tube_axis * along;
    m.exit_lateral_velocity = lateral.length();

    per_seed_.push_back(m);
}

static double mean(const std::vector<double>& v) {
    if (v.empty()) return 0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

static double stddev(const std::vector<double>& v, double m) {
    if (v.size() < 2) return 0;
    double sum = 0;
    for (double x : v) sum += (x - m) * (x - m);
    return std::sqrt(sum / (v.size() - 1));
}

void MetricsCollector::compute_summary() {
    summary_ = {};
    summary_.total_seeds = (int)per_seed_.size();
    summary_.exited_seeds = (int)per_seed_.size();

    if (per_seed_.empty()) return;

    std::vector<double> times, hits, speeds, laterals, impulses, bounces;
    for (auto& m : per_seed_) {
        times.push_back(m.time_to_exit);
        hits.push_back(m.wall_hits);
        speeds.push_back(m.exit_speed);
        laterals.push_back(m.exit_lateral_velocity);
        impulses.push_back(m.cumulative_impulse);
        bounces.push_back(m.max_bounce_energy);
    }

    summary_.mean_time_to_exit = mean(times);
    summary_.sd_time_to_exit = stddev(times, summary_.mean_time_to_exit);
    summary_.mean_wall_hits = mean(hits);
    summary_.sd_wall_hits = stddev(hits, summary_.mean_wall_hits);
    summary_.mean_exit_speed = mean(speeds);
    summary_.sd_exit_speed = stddev(speeds, summary_.mean_exit_speed);
    summary_.mean_exit_lateral_vel = mean(laterals);
    summary_.sd_exit_lateral_vel = stddev(laterals, summary_.mean_exit_lateral_vel);
    summary_.mean_cumulative_impulse = mean(impulses);
    summary_.mean_max_bounce_energy = mean(bounces);

    // Spacing risk proxy: SD(time_to_exit) + SD(exit_lateral_vel) scaled
    // Scale lateral vel SD to same order as time SD
    double time_scale = (summary_.mean_time_to_exit > 0) ? 1.0 / summary_.mean_time_to_exit : 1.0;
    double lat_scale = (summary_.mean_exit_speed > 0) ? 1.0 / summary_.mean_exit_speed : 1.0;
    summary_.spacing_risk_proxy = summary_.sd_time_to_exit * time_scale
                                 + summary_.sd_exit_lateral_vel * lat_scale;
}

bool MetricsCollector::write_csv(const std::string& path) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;

    f << "id,time_to_exit,wall_hits,cumulative_impulse,max_bounce_energy,"
         "exit_vx,exit_vy,exit_vz,exit_speed,exit_lateral_velocity\n";
    f << std::fixed << std::setprecision(6);

    for (auto& m : per_seed_) {
        f << m.id << ","
          << m.time_to_exit << ","
          << m.wall_hits << ","
          << m.cumulative_impulse << ","
          << m.max_bounce_energy << ","
          << m.exit_velocity.x << ","
          << m.exit_velocity.y << ","
          << m.exit_velocity.z << ","
          << m.exit_speed << ","
          << m.exit_lateral_velocity << "\n";
    }
    return true;
}

bool MetricsCollector::write_summary(const std::string& path) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;

    f << std::fixed << std::setprecision(6);
    f << "metric,value\n";
    f << "total_seeds," << summary_.total_seeds << "\n";
    f << "exited_seeds," << summary_.exited_seeds << "\n";
    f << "mean_time_to_exit," << summary_.mean_time_to_exit << "\n";
    f << "sd_time_to_exit," << summary_.sd_time_to_exit << "\n";
    f << "mean_wall_hits," << summary_.mean_wall_hits << "\n";
    f << "sd_wall_hits," << summary_.sd_wall_hits << "\n";
    f << "mean_exit_speed," << summary_.mean_exit_speed << "\n";
    f << "sd_exit_speed," << summary_.sd_exit_speed << "\n";
    f << "mean_exit_lateral_vel," << summary_.mean_exit_lateral_vel << "\n";
    f << "sd_exit_lateral_vel," << summary_.sd_exit_lateral_vel << "\n";
    f << "mean_cumulative_impulse," << summary_.mean_cumulative_impulse << "\n";
    f << "mean_max_bounce_energy," << summary_.mean_max_bounce_energy << "\n";
    f << "spacing_risk_proxy," << summary_.spacing_risk_proxy << "\n";
    return true;
}

bool MetricsCollector::write_vtk_trajectories(const std::string& path,
                                                const std::vector<std::vector<vec3>>& trajectories,
                                                const std::vector<std::vector<double>>& times) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;

    // Count total points and lines
    int total_points = 0;
    for (auto& t : trajectories) total_points += (int)t.size();

    f << "# vtk DataFile Version 3.0\n";
    f << "DEM seed trajectories\n";
    f << "ASCII\n";
    f << "DATASET POLYDATA\n";
    f << "POINTS " << total_points << " double\n";

    for (auto& traj : trajectories) {
        for (auto& p : traj) {
            f << p.x << " " << p.y << " " << p.z << "\n";
        }
    }

    // Lines
    int n_lines = (int)trajectories.size();
    int total_indices = total_points + n_lines;
    f << "LINES " << n_lines << " " << total_indices << "\n";

    int offset = 0;
    for (auto& traj : trajectories) {
        f << traj.size();
        for (size_t i = 0; i < traj.size(); i++) {
            f << " " << (offset + i);
        }
        f << "\n";
        offset += (int)traj.size();
    }

    // Point data: time
    f << "POINT_DATA " << total_points << "\n";
    f << "SCALARS time double\n";
    f << "LOOKUP_TABLE default\n";
    for (auto& t : times) {
        for (double v : t) f << v << "\n";
    }

    return true;
}

bool MetricsCollector::write_trajectory_csv(const std::string& path,
                                              const std::vector<std::vector<vec3>>& trajectories,
                                              const std::vector<std::vector<double>>& times,
                                              const std::vector<uint32_t>& ids) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;

    f << "particle_id,time,x,y,z\n";
    f << std::fixed << std::setprecision(6);

    for (size_t i = 0; i < trajectories.size(); i++) {
        uint32_t id = (i < ids.size()) ? ids[i] : (uint32_t)i;
        for (size_t j = 0; j < trajectories[i].size(); j++) {
            f << id << ","
              << times[i][j] << ","
              << trajectories[i][j].x << ","
              << trajectories[i][j].y << ","
              << trajectories[i][j].z << "\n";
        }
    }
    return true;
}

// ─── Trajectory Recorder ─────────────────────────────────────────────────────

void TrajectoryRecorder::init(int max_track) {
    max_particles_to_track = max_track;
    tracked_ids_.clear();
    positions_.clear();
    times_.clear();
    last_sample_time_ = -1e30;
}

void TrajectoryRecorder::sample(const ParticleSystem& psys, double current_time) {
    if (current_time - last_sample_time_ < sample_interval) return;
    last_sample_time_ = current_time;

    // Track new particles
    for (auto& p : psys.particles) {
        if (!p.active) continue;
        if ((int)tracked_ids_.size() >= max_particles_to_track) break;

        bool found = false;
        for (auto id : tracked_ids_) {
            if (id == p.id) { found = true; break; }
        }
        if (!found) {
            tracked_ids_.push_back(p.id);
            positions_.push_back({});
            times_.push_back({});
        }
    }

    // Record positions
    for (size_t i = 0; i < tracked_ids_.size(); i++) {
        for (auto& p : psys.particles) {
            if (p.id == tracked_ids_[i] && p.active) {
                positions_[i].push_back(p.pos);
                times_[i].push_back(current_time);
                break;
            }
        }
    }
}

} // namespace dem
