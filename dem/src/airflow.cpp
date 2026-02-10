// Airflow field: LBM velocity grid reader and interpolation

#include "airflow.hpp"
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace dem {

bool AirflowField::load_binary(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    // Header: nx, ny, nz (int32), origin x,y,z (double), spacing (double)
    int32_t dims[3];
    file.read(reinterpret_cast<char*>(dims), sizeof(dims));
    nx = dims[0]; ny = dims[1]; nz = dims[2];

    double orig[3];
    file.read(reinterpret_cast<char*>(orig), sizeof(orig));
    origin = {orig[0], orig[1], orig[2]};

    file.read(reinterpret_cast<char*>(&spacing), sizeof(spacing));

    if (!file || nx <= 0 || ny <= 0 || nz <= 0) return false;

    size_t total = (size_t)nx * ny * nz;
    ux.resize(total);
    uy.resize(total);
    uz.resize(total);

    file.read(reinterpret_cast<char*>(ux.data()), total * sizeof(float));
    file.read(reinterpret_cast<char*>(uy.data()), total * sizeof(float));
    file.read(reinterpret_cast<char*>(uz.data()), total * sizeof(float));

    // Optional density
    density.resize(total, (float)rho_fluid);
    file.read(reinterpret_cast<char*>(density.data()), total * sizeof(float));
    // It's OK if density read fails - we already filled with default

    return file.good() || file.eof();
}

bool AirflowField::load_npy(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    // Minimal NumPy .npy format parser
    char magic[6];
    file.read(magic, 6);
    if (std::memcmp(magic, "\x93NUMPY", 6) != 0) return false;

    uint8_t major, minor;
    file.read(reinterpret_cast<char*>(&major), 1);
    file.read(reinterpret_cast<char*>(&minor), 1);

    uint16_t header_len;
    if (major == 1) {
        file.read(reinterpret_cast<char*>(&header_len), 2);
    } else {
        uint32_t hl;
        file.read(reinterpret_cast<char*>(&hl), 4);
        header_len = (uint16_t)hl;
    }

    std::string header(header_len, '\0');
    file.read(&header[0], header_len);

    // Parse shape from header: 'shape': (3, nz, ny, nx)
    // Simple parsing
    auto shape_pos = header.find("'shape'");
    if (shape_pos == std::string::npos) shape_pos = header.find("\"shape\"");
    if (shape_pos == std::string::npos) return false;

    auto paren_start = header.find('(', shape_pos);
    auto paren_end = header.find(')', shape_pos);
    if (paren_start == std::string::npos || paren_end == std::string::npos) return false;

    std::string shape_str = header.substr(paren_start + 1, paren_end - paren_start - 1);
    int dims[4] = {0};
    int dim_idx = 0;
    std::string num;
    for (char c : shape_str) {
        if (c >= '0' && c <= '9') {
            num += c;
        } else if (!num.empty()) {
            if (dim_idx < 4) dims[dim_idx++] = std::stoi(num);
            num.clear();
        }
    }
    if (!num.empty() && dim_idx < 4) dims[dim_idx++] = std::stoi(num);

    if (dim_idx < 4 || dims[0] != 3) return false;
    nz = dims[1]; ny = dims[2]; nx = dims[3];

    // Check dtype
    bool is_float32 = header.find("<f4") != std::string::npos ||
                       header.find("float32") != std::string::npos;
    bool is_float64 = header.find("<f8") != std::string::npos ||
                       header.find("float64") != std::string::npos;

    size_t total = (size_t)nx * ny * nz;
    ux.resize(total);
    uy.resize(total);
    uz.resize(total);

    if (is_float32) {
        file.read(reinterpret_cast<char*>(ux.data()), total * sizeof(float));
        file.read(reinterpret_cast<char*>(uy.data()), total * sizeof(float));
        file.read(reinterpret_cast<char*>(uz.data()), total * sizeof(float));
    } else if (is_float64) {
        std::vector<double> buf(total);
        file.read(reinterpret_cast<char*>(buf.data()), total * sizeof(double));
        for (size_t i = 0; i < total; i++) ux[i] = (float)buf[i];
        file.read(reinterpret_cast<char*>(buf.data()), total * sizeof(double));
        for (size_t i = 0; i < total; i++) uy[i] = (float)buf[i];
        file.read(reinterpret_cast<char*>(buf.data()), total * sizeof(double));
        for (size_t i = 0; i < total; i++) uz[i] = (float)buf[i];
    } else {
        return false;
    }

    density.resize(total, (float)rho_fluid);
    return true;
}

void AirflowField::world_to_grid(const vec3& pos, double& gx, double& gy, double& gz) const {
    gx = (pos.x - origin.x) / spacing;
    gy = (pos.y - origin.y) / spacing;
    gz = (pos.z - origin.z) / spacing;
}

float AirflowField::trilinear(const std::vector<float>& data, double gx, double gy, double gz) const {
    int ix = (int)std::floor(gx);
    int iy = (int)std::floor(gy);
    int iz = (int)std::floor(gz);

    double fx = gx - ix;
    double fy = gy - iy;
    double fz = gz - iz;

    // Clamp to valid range
    auto clamp_idx = [](int i, int max_val) -> int {
        return std::max(0, std::min(i, max_val - 1));
    };

    int x0 = clamp_idx(ix, nx), x1 = clamp_idx(ix + 1, nx);
    int y0 = clamp_idx(iy, ny), y1 = clamp_idx(iy + 1, ny);
    int z0 = clamp_idx(iz, nz), z1 = clamp_idx(iz + 1, nz);

    auto idx = [&](int x, int y, int z) -> size_t {
        return (size_t)x + nx * ((size_t)y + ny * (size_t)z);
    };

    float c000 = data[idx(x0,y0,z0)], c100 = data[idx(x1,y0,z0)];
    float c010 = data[idx(x0,y1,z0)], c110 = data[idx(x1,y1,z0)];
    float c001 = data[idx(x0,y0,z1)], c101 = data[idx(x1,y0,z1)];
    float c011 = data[idx(x0,y1,z1)], c111 = data[idx(x1,y1,z1)];

    float c00 = c000 * (1-fx) + c100 * fx;
    float c01 = c001 * (1-fx) + c101 * fx;
    float c10 = c010 * (1-fx) + c110 * fx;
    float c11 = c011 * (1-fx) + c111 * fx;

    float c0 = c00 * (1-fy) + c10 * fy;
    float c1 = c01 * (1-fy) + c11 * fy;

    return c0 * (1-fz) + c1 * fz;
}

vec3 AirflowField::interpolate_velocity(const vec3& pos) const {
    if (!is_loaded()) return {0,0,0};

    double gx, gy, gz;
    world_to_grid(pos, gx, gy, gz);

    // Check bounds (with a little margin)
    if (gx < -0.5 || gx > nx - 0.5 || gy < -0.5 || gy > ny - 0.5 || gz < -0.5 || gz > nz - 0.5)
        return {0,0,0};

    return {
        (double)trilinear(ux, gx, gy, gz),
        (double)trilinear(uy, gx, gy, gz),
        (double)trilinear(uz, gx, gy, gz)
    };
}

double AirflowField::interpolate_density(const vec3& pos) const {
    if (!is_loaded() || density.empty()) return rho_fluid;
    double gx, gy, gz;
    world_to_grid(pos, gx, gy, gz);
    return (double)trilinear(density, gx, gy, gz);
}

vec3 AirflowField::compute_drag(const vec3& fluid_vel, const vec3& particle_vel,
                                  double rho_fluid, double Cd, double particle_radius) {
    vec3 v_rel = fluid_vel - particle_vel;
    double v_rel_mag = v_rel.length();
    if (v_rel_mag < 1e-12) return {0,0,0};

    double A = M_PI * particle_radius * particle_radius;
    double Fd_mag = 0.5 * rho_fluid * Cd * A * v_rel_mag * v_rel_mag;
    return v_rel.normalized() * Fd_mag;
}

double AirflowField::schiller_naumann_Cd(double Re) {
    if (Re < 1e-6) return 0.0;
    if (Re <= 1000.0) {
        return (24.0 / Re) * (1.0 + 0.15 * std::pow(Re, 0.687));
    }
    return 0.44;  // Newton regime
}

} // namespace dem
