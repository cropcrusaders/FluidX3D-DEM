// Seed injector implementation

#include "injector.hpp"
#include <cmath>

namespace dem {

void Injector::init(const InjectorConfig& cfg, RNG& rng) {
    config = cfg;
    total_injected_ = 0;
    next_inject_time_ = 0.0;
}

int Injector::inject(ParticleSystem& psys, double current_time, double dt, RNG& rng) {
    if (done()) return 0;

    int count = 0;
    double interval = 1.0 / config.seeds_per_second;
    if (config.seed_spacing_time > 0) {
        interval = std::max(interval, config.seed_spacing_time);
    }

    while (current_time >= next_inject_time_ && !done()) {
        // Generate position within aperture
        vec3 pos = config.position;

        if (config.aperture == ApertureShape::CIRCULAR) {
            // Uniform random in circle
            double r = config.aperture_radius * std::sqrt(rng.uniform());
            double theta = rng.uniform(0, 2.0 * M_PI);

            // Build local frame perpendicular to injection direction
            vec3 dir = config.direction.normalized();
            vec3 up = (std::abs(dir.y) < 0.9) ? vec3(0,1,0) : vec3(1,0,0);
            vec3 right = cross(dir, up).normalized();
            up = cross(right, dir).normalized();

            pos = pos + right * (r * std::cos(theta)) + up * (r * std::sin(theta));
        } else {
            // Rectangular aperture
            double ox = rng.uniform(-0.5, 0.5) * config.aperture_width;
            double oy = rng.uniform(-0.5, 0.5) * config.aperture_height;

            vec3 dir = config.direction.normalized();
            vec3 up = (std::abs(dir.y) < 0.9) ? vec3(0,1,0) : vec3(1,0,0);
            vec3 right = cross(dir, up).normalized();
            up = cross(right, dir).normalized();

            pos = pos + right * ox + up * oy;
        }

        // Generate velocity
        double speed = config.initial_speed;
        if (config.speed_stddev > 0) {
            speed += rng.normal(0, config.speed_stddev);
        }
        vec3 vel = config.direction.normalized() * speed;

        // Add lateral scatter
        if (config.lateral_speed_stddev > 0) {
            vec3 dir = config.direction.normalized();
            vec3 up = (std::abs(dir.y) < 0.9) ? vec3(0,1,0) : vec3(1,0,0);
            vec3 right = cross(dir, up).normalized();
            up = cross(right, dir).normalized();

            vel = vel + right * rng.normal(0, config.lateral_speed_stddev)
                      + up * rng.normal(0, config.lateral_speed_stddev);
        }

        psys.add_particle(pos, vel, config.seed_type_idx, current_time);
        total_injected_++;
        count++;
        next_inject_time_ += interval;
    }

    return count;
}

} // namespace dem
