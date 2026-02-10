#pragma once
// Seed injector: generates particles at the tube inlet with configurable rate/distribution

#include "dem_math.hpp"
#include "particle.hpp"
#include <string>

namespace dem {

enum class ApertureShape { CIRCULAR, RECTANGULAR };

struct InjectorConfig {
    // Position
    vec3 position = {0, 0, 0.3};  // inlet center
    vec3 direction = {0, 0, -1};  // injection direction (normalized)

    // Rate
    double seeds_per_second = 100.0;
    int max_seeds = 5000;

    // Aperture
    ApertureShape aperture = ApertureShape::CIRCULAR;
    double aperture_radius = 0.02;   // m
    double aperture_width = 0.02;    // m (rectangular)
    double aperture_height = 0.02;   // m (rectangular)

    // Initial velocity
    double initial_speed = 0.0;          // m/s in injection direction
    double speed_stddev = 0.0;           // m/s
    double lateral_speed_stddev = 0.01;  // m/s lateral scatter

    // Seed type
    int seed_type_idx = 0;

    // Spacing (singulation)
    double seed_spacing_time = 0.0;  // extra time between seeds (0 = use rate only)
};

class Injector {
public:
    InjectorConfig config;

    void init(const InjectorConfig& cfg, RNG& rng);

    // Call each step; returns number of seeds injected this step
    int inject(ParticleSystem& psys, double current_time, double dt, RNG& rng);

    int total_injected() const { return total_injected_; }
    bool done() const { return total_injected_ >= config.max_seeds; }

private:
    double next_inject_time_ = 0.0;
    int total_injected_ = 0;
};

// Exit plane: detects when particles cross and records metrics
struct ExitPlane {
    vec3 point = {0, 0, 0};       // point on plane
    vec3 normal = {0, 0, 1};      // outward normal (particles exit when crossing in -normal dir)

    // Check if particle has crossed the exit plane
    bool has_crossed(const vec3& pos) const {
        return dot(pos - point, normal) < 0;
    }
};

} // namespace dem
