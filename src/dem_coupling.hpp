#pragma once

#include "defines.hpp"

#ifdef DEM_COUPLING

#include "../dem/include/simulator.hpp"
#include <mutex>
#include <vector>

struct DemParticleVis {
	float x, y, z; // position in LBM centered coordinates (lattice units)
	float radius;   // radius in lattice units
	int color;
};

class DemCoupling {
public:
	dem::Simulator sim;

	// Unit conversion
	float spacing = 0.001f;     // meters per LBM cell = units.si_x(1.0f)
	float u_conversion = 1.0f;  // LBM velocity to SI = units.si_u(1.0f)
	float dt_lbm_si = 1e-4f;    // LBM timestep in seconds = units.si_t(1)
	uint Nx=0u, Ny=0u, Nz=0u;

	// Thread-safe visualization buffer
	std::mutex vis_mutex;
	std::vector<DemParticleVis> vis_particles;

	// Velocity buffer for coupling (SI units)
	std::vector<float> ux_buf, uy_buf, uz_buf;

	// Initialize coupling: set up DEM airflow grid to match LBM grid
	void init(uint Nx, uint Ny, uint Nz, float spacing, float u_conversion, float dt_lbm_si);

	// Transfer LBM velocity field to DEM and step DEM forward by one LBM timestep
	// Reads velocity from LBM (GPU->CPU), converts to SI, feeds to DEM airflow,
	// then runs DEM substeps to advance by dt_lbm_si seconds
	void step(class LBM& lbm);

	// Update visualization buffer from current DEM particle positions
	void update_vis();

	// Draw DEM particles onto screen (called from graphics thread)
	void draw_particles();
};

extern DemCoupling* g_dem_coupling;

#endif // DEM_COUPLING
