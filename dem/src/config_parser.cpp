// JSON config parser implementation

#include "config_parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>

namespace dem {

// ─── Minimal JSON parser ─────────────────────────────────────────────────────

static void skip_ws(const std::string& s, size_t& i) {
    while (i < s.size() && std::isspace(s[i])) i++;
}

static JsonValue parse_value(const std::string& s, size_t& i);

static std::string parse_string(const std::string& s, size_t& i) {
    if (s[i] != '"') return "";
    i++;
    std::string result;
    while (i < s.size() && s[i] != '"') {
        if (s[i] == '\\') {
            i++;
            if (i < s.size()) {
                switch (s[i]) {
                    case '"': result += '"'; break;
                    case '\\': result += '\\'; break;
                    case '/': result += '/'; break;
                    case 'n': result += '\n'; break;
                    case 't': result += '\t'; break;
                    case 'r': result += '\r'; break;
                    default: result += s[i]; break;
                }
            }
        } else {
            result += s[i];
        }
        i++;
    }
    if (i < s.size()) i++;  // skip closing quote
    return result;
}

static double parse_number(const std::string& s, size_t& i) {
    size_t start = i;
    if (s[i] == '-') i++;
    while (i < s.size() && (std::isdigit(s[i]) || s[i] == '.' || s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-')) {
        if ((s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-') && i == start) break;
        i++;
    }
    return std::stod(s.substr(start, i - start));
}

static JsonValue parse_value(const std::string& s, size_t& i) {
    skip_ws(s, i);
    if (i >= s.size()) return {};

    if (s[i] == '"') {
        return JsonValue(parse_string(s, i));
    }
    if (s[i] == '{') {
        i++;
        JsonValue obj;
        skip_ws(s, i);
        if (i < s.size() && s[i] == '}') { i++; obj.set("__dummy__", JsonValue()); return obj; }
        while (i < s.size()) {
            skip_ws(s, i);
            if (s[i] == '}') { i++; break; }
            if (s[i] == ',') { i++; continue; }
            std::string key = parse_string(s, i);
            skip_ws(s, i);
            if (i < s.size() && s[i] == ':') i++;
            JsonValue val = parse_value(s, i);
            obj.set(key, val);
            skip_ws(s, i);
        }
        return obj;
    }
    if (s[i] == '[') {
        i++;
        JsonValue arr;
        skip_ws(s, i);
        if (i < s.size() && s[i] == ']') { i++; return arr; }
        while (i < s.size()) {
            skip_ws(s, i);
            if (s[i] == ']') { i++; break; }
            if (s[i] == ',') { i++; continue; }
            arr.push_back(parse_value(s, i));
            skip_ws(s, i);
        }
        return arr;
    }
    if (s[i] == 't') { i += 4; return JsonValue(true); }
    if (s[i] == 'f') { i += 5; return JsonValue(false); }
    if (s[i] == 'n') { i += 4; return JsonValue(); }

    // Number
    return JsonValue(parse_number(s, i));
}

JsonValue JsonValue::parse(const std::string& json) {
    // Strip comments (// style)
    std::string clean;
    bool in_string = false;
    for (size_t i = 0; i < json.size(); i++) {
        if (json[i] == '"' && (i == 0 || json[i-1] != '\\')) {
            in_string = !in_string;
        }
        if (!in_string && i + 1 < json.size() && json[i] == '/' && json[i+1] == '/') {
            while (i < json.size() && json[i] != '\n') i++;
            continue;
        }
        clean += json[i];
    }

    size_t pos = 0;
    return parse_value(clean, pos);
}

// ─── Config loading ──────────────────────────────────────────────────────────

static vec3 parse_vec3(const JsonValue& v, vec3 def = {0,0,0}) {
    if (v.is_array() && v.size() >= 3) {
        return {v[0].as_number(), v[1].as_number(), v[2].as_number()};
    }
    return def;
}

bool load_config(const std::string& path, Simulator& sim) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Error: Cannot open config file: " << path << "\n";
        return false;
    }
    std::stringstream ss;
    ss << f.rdbuf();
    return parse_config(ss.str(), sim);
}

bool parse_config(const std::string& json_str, Simulator& sim) {
    auto root = JsonValue::parse(json_str);
    if (!root.is_object()) {
        std::cerr << "Error: Config root must be a JSON object\n";
        return false;
    }

    auto& cfg = sim.config;

    // Simulation settings
    if (root.has("simulation")) {
        auto& s = root["simulation"];
        cfg.dt = s["dt"].as_number(cfg.dt);
        cfg.max_time = s["max_time"].as_number(cfg.max_time);
        cfg.output_interval = s["output_interval"].as_int(cfg.output_interval);
        cfg.gravity = parse_vec3(s["gravity"], cfg.gravity);
        cfg.auto_dt = s["auto_dt"].as_bool(cfg.auto_dt);
        cfg.dt_safety = s["dt_safety"].as_number(cfg.dt_safety);
        cfg.rng_seed = (uint64_t)s["rng_seed"].as_number((double)cfg.rng_seed);
        cfg.output_dir = s["output_dir"].as_string(cfg.output_dir);
        cfg.save_trajectories = s["save_trajectories"].as_bool(cfg.save_trajectories);
        cfg.max_trajectory_particles = s["max_trajectory_particles"].as_int(cfg.max_trajectory_particles);
        cfg.trajectory_sample_interval = s["trajectory_sample_interval"].as_number(cfg.trajectory_sample_interval);
        cfg.tube_axis = parse_vec3(s["tube_axis"], cfg.tube_axis);
    }

    // Seed types
    if (root.has("seed_types")) {
        auto& types = root["seed_types"];
        for (size_t i = 0; i < types.size(); i++) {
            auto& t = types[i];
            SeedType st;
            st.name = t["name"].as_string("default");
            std::string shape = t["shape"].as_string("sphere");
            if (shape == "ellipsoid" || shape == "multi_sphere_ellipsoid") {
                st.shape = SeedShape::MULTI_SPHERE_ELLIPSOID;
            }
            st.radius = t["radius"].as_number(st.radius);
            st.mass = t["mass"].as_number(st.mass);
            st.a = t["a"].as_number(st.radius);
            st.b = t["b"].as_number(st.radius * 0.75);
            st.c = t["c"].as_number(st.radius * 0.5);
            sim.particles.seed_types.push_back(st);
        }
    }
    // Default seed type if none specified
    if (sim.particles.seed_types.empty()) {
        SeedType st;
        st.name = "default";
        sim.particles.seed_types.push_back(st);
    }

    // Contact properties
    if (root.has("contact")) {
        auto& c = root["contact"];
        if (c.has("seed_seed")) {
            auto& ss = c["seed_seed"];
            sim.particles.default_seed_seed.kn = ss["kn"].as_number(1e4);
            sim.particles.default_seed_seed.kt = ss["kt"].as_number(2500);
            sim.particles.default_seed_seed.restitution = ss["restitution"].as_number(0.5);
            sim.particles.default_seed_seed.friction_s = ss["friction"].as_number(0.5);
            sim.particles.default_seed_seed.friction_r = ss["rolling_friction"].as_number(0.01);
        }
        if (c.has("seed_wall")) {
            auto& sw = c["seed_wall"];
            sim.particles.default_seed_wall.kn = sw["kn"].as_number(1e4);
            sim.particles.default_seed_wall.kt = sw["kt"].as_number(2500);
            sim.particles.default_seed_wall.restitution = sw["restitution"].as_number(0.4);
            sim.particles.default_seed_wall.friction_s = sw["friction"].as_number(0.4);
            sim.particles.default_seed_wall.friction_r = sw["rolling_friction"].as_number(0.01);
        }
    } else {
        // Defaults
        sim.particles.default_seed_seed = {1e4, 2500, 0.5, 0.5, 0.01};
        sim.particles.default_seed_wall = {1e4, 2500, 0.4, 0.4, 0.01};
    }

    // Geometry
    if (root.has("geometry")) {
        auto& g = root["geometry"];

        // Mesh files
        if (g.has("mesh")) {
            std::string mesh_path = g["mesh"].as_string();
            if (!mesh_path.empty()) {
                if (!sim.geometry.load_stl(mesh_path)) {
                    std::cerr << "Warning: Failed to load mesh: " << mesh_path << "\n";
                }
            }
        }

        // Primitives
        if (g.has("primitives")) {
            auto& prims = g["primitives"];
            for (size_t i = 0; i < prims.size(); i++) {
                auto& pr = prims[i];
                AnalyticPrimitive prim;
                std::string type = pr["type"].as_string("cylinder");

                if (type == "cylinder") prim.type = PrimitiveType::CYLINDER;
                else if (type == "cone" || type == "frustum") prim.type = PrimitiveType::CONE;
                else if (type == "bent_cylinder" || type == "curve") prim.type = PrimitiveType::BENT_CYLINDER;
                else if (type == "spiral" || type == "spiral_insert") prim.type = PrimitiveType::SPIRAL_INSERT;
                else if (type == "plane") prim.type = PrimitiveType::PLANE;

                prim.center = parse_vec3(pr["center"], prim.center);
                prim.axis = parse_vec3(pr["axis"], prim.axis);
                prim.length = pr["length"].as_number(prim.length);
                prim.radius = pr["radius"].as_number(prim.radius);
                prim.radius_top = pr["radius_top"].as_number(prim.radius);
                prim.radius_bottom = pr["radius_bottom"].as_number(prim.radius);
                prim.bend_radius = pr["bend_radius"].as_number(prim.bend_radius);
                prim.bend_angle = pr["bend_angle"].as_number(prim.bend_angle);
                prim.spiral_pitch = pr["spiral_pitch"].as_number(prim.spiral_pitch);
                prim.spiral_height = pr["spiral_height"].as_number(prim.spiral_height);
                prim.spiral_turns = pr["spiral_turns"].as_int(prim.spiral_turns);
                prim.plane_normal = parse_vec3(pr["plane_normal"], prim.plane_normal);
                prim.plane_offset = pr["plane_offset"].as_number(prim.plane_offset);

                sim.geometry.add_primitive(prim);
            }
        }
    }

    // Injector
    if (root.has("injector")) {
        auto& inj = root["injector"];
        auto& ic = sim.injector.config;
        ic.position = parse_vec3(inj["position"], ic.position);
        ic.direction = parse_vec3(inj["direction"], ic.direction);
        ic.seeds_per_second = inj["seeds_per_second"].as_number(ic.seeds_per_second);
        ic.max_seeds = inj["max_seeds"].as_int(ic.max_seeds);

        std::string ap = inj["aperture_shape"].as_string("circular");
        if (ap == "rectangular") ic.aperture = ApertureShape::RECTANGULAR;

        ic.aperture_radius = inj["aperture_radius"].as_number(ic.aperture_radius);
        ic.aperture_width = inj["aperture_width"].as_number(ic.aperture_width);
        ic.aperture_height = inj["aperture_height"].as_number(ic.aperture_height);
        ic.initial_speed = inj["initial_speed"].as_number(ic.initial_speed);
        ic.speed_stddev = inj["speed_stddev"].as_number(ic.speed_stddev);
        ic.lateral_speed_stddev = inj["lateral_speed_stddev"].as_number(ic.lateral_speed_stddev);
        ic.seed_type_idx = inj["seed_type_idx"].as_int(ic.seed_type_idx);
        ic.seed_spacing_time = inj["seed_spacing_time"].as_number(ic.seed_spacing_time);
    }

    // Exit plane
    if (root.has("exit_plane")) {
        auto& ep = root["exit_plane"];
        sim.exit_plane.point = parse_vec3(ep["point"], sim.exit_plane.point);
        sim.exit_plane.normal = parse_vec3(ep["normal"], sim.exit_plane.normal);
    }

    // Airflow
    if (root.has("airflow")) {
        auto& af = root["airflow"];
        cfg.enable_airflow = af["enabled"].as_bool(false);
        cfg.airflow_path = af["path"].as_string("");
        cfg.airflow_format = af["format"].as_string("binary");
        cfg.drag_Cd = af["drag_Cd"].as_number(cfg.drag_Cd);
        cfg.use_schiller_naumann = af["schiller_naumann"].as_bool(false);
        cfg.air_density = af["air_density"].as_number(cfg.air_density);

        if (af.has("grid")) {
            sim.airflow.origin = parse_vec3(af["grid"]["origin"], sim.airflow.origin);
            sim.airflow.spacing = af["grid"]["spacing"].as_number(sim.airflow.spacing);
        }
    }

    return true;
}

} // namespace dem
