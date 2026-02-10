#pragma once
// LBM velocity field reader for one-way air-seed coupling
// Reads a 3D uniform grid of velocity data exported from FluidX3D

#include "dem_math.hpp"
#include <string>
#include <vector>

namespace dem {

class AirflowField {
public:
    // Grid parameters
    vec3 origin = {0,0,0};
    double spacing = 0.001;  // grid cell size (m)
    int nx = 0, ny = 0, nz = 0;

    // Fluid properties
    double rho_fluid = 1.225;  // air density (kg/m^3)

    // Velocity data (linearized: index = x + nx*(y + ny*z))
    std::vector<float> ux, uy, uz;
    std::vector<float> density;  // optional

    // Load from raw binary file
    // Format: header (3 ints: nx,ny,nz + 3 doubles: ox,oy,oz + 1 double: spacing)
    // Then nx*ny*nz floats for ux, then uy, then uz
    bool load_binary(const std::string& path);

    // Load from NumPy .npy file (4D array: [3, nz, ny, nx])
    bool load_npy(const std::string& path);

    // Trilinear interpolation of velocity at world position
    vec3 interpolate_velocity(const vec3& pos) const;

    // Interpolate density at world position
    double interpolate_density(const vec3& pos) const;

    bool is_loaded() const { return nx > 0 && ny > 0 && nz > 0 && !ux.empty(); }

    // Compute drag force on a sphere
    // Fd = 0.5 * rho * Cd * A * |v_rel| * v_rel
    static vec3 compute_drag(const vec3& fluid_vel, const vec3& particle_vel,
                             double rho_fluid, double Cd, double particle_radius);

    // Schiller-Naumann drag coefficient
    static double schiller_naumann_Cd(double Re);

private:
    // Convert world position to grid coordinates
    void world_to_grid(const vec3& pos, double& gx, double& gy, double& gz) const;

    // Trilinear interpolation helper
    float trilinear(const std::vector<float>& data, double gx, double gy, double gz) const;
};

} // namespace dem
