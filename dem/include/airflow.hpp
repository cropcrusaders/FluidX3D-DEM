#pragma once
// LBM velocity field reader for air-seed coupling (one-way and two-way)
// Reads a 3D uniform grid of velocity data exported from FluidX3D,
// or accepts live velocity updates for two-way coupled simulations.

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

    // Initialize grid dimensions without loading data (for two-way coupling)
    void init_grid(int nx_, int ny_, int nz_, const vec3& origin_, double spacing_);

    // Update velocity field from external float arrays (for live two-way coupling)
    // Arrays must have nx*ny*nz elements, linearized as index = x + nx*(y + ny*z)
    // Velocity values must be in SI units (m/s)
    void update_from_arrays(const float* ux_data, const float* uy_data, const float* uz_data);

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

    // ─── Two-way coupling: force feedback ─────────────────────────────
    // Per-particle drag force record for feeding back to LBM
    struct ParticleDrag {
        vec3 pos;       // particle position (world space)
        vec3 drag;      // drag force on particle (N) -- reaction = -drag goes to fluid
        double radius;  // particle radius (m)
    };

    // Distribute a point force onto the nearest grid cells using trilinear weighting
    // Writes force contributions into fx_out, fy_out, fz_out arrays (must be nx*ny*nz)
    // The force applied to the fluid is the reaction: -particle_drag (Newton's 3rd law)
    void distribute_force_to_grid(const ParticleDrag& pd,
                                  float* fx_out, float* fy_out, float* fz_out) const;

    // Distribute forces from all particles in the list
    void distribute_forces_to_grid(const std::vector<ParticleDrag>& drags,
                                   float* fx_out, float* fy_out, float* fz_out) const;

private:
    // Convert world position to grid coordinates
    void world_to_grid(const vec3& pos, double& gx, double& gy, double& gz) const;

    // Trilinear interpolation helper
    float trilinear(const std::vector<float>& data, double gx, double gy, double gz) const;
};

} // namespace dem
