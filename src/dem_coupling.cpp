#include "dem_coupling.hpp"

#ifdef DEM_COUPLING

#include "lbm.hpp"
#include "graphics.hpp"
#include <cstring>
#include <iostream>
#include <cmath>

DemCoupling* g_dem_coupling = nullptr;

void DemCoupling::init(uint Nx_, uint Ny_, uint Nz_, float spacing_, float u_conversion_, float dt_lbm_si_, float f_conversion_) {
	Nx = Nx_;
	Ny = Ny_;
	Nz = Nz_;
	spacing = spacing_;
	u_conversion = u_conversion_;
	dt_lbm_si = dt_lbm_si_;
	f_conversion = f_conversion_;

	// Allocate velocity buffers
	const ulong N = (ulong)Nx * (ulong)Ny * (ulong)Nz;
	ux_buf.resize(N, 0.0f);
	uy_buf.resize(N, 0.0f);
	uz_buf.resize(N, 0.0f);

	// Allocate force feedback buffers for two-way coupling
	fx_buf.resize(N, 0.0f);
	fy_buf.resize(N, 0.0f);
	fz_buf.resize(N, 0.0f);

	// Initialize DEM airflow grid to match LBM grid
	// LBM cell (ix,iy,iz) has centered position: (ix-Nx/2+0.5, iy-Ny/2+0.5, iz-Nz/2+0.5) * spacing
	// DEM grid cell (ix,iy,iz) has position: origin + (ix,iy,iz) * spacing
	// So origin = (-Nx/2+0.5, -Ny/2+0.5, -Nz/2+0.5) * spacing
	dem::vec3 origin(
		(-(double)Nx/2.0 + 0.5) * (double)spacing,
		(-(double)Ny/2.0 + 0.5) * (double)spacing,
		(-(double)Nz/2.0 + 0.5) * (double)spacing
	);
	sim.airflow.init_grid((int)Nx, (int)Ny, (int)Nz, origin, (double)spacing);
	sim.airflow.rho_fluid = sim.config.air_density;

	// Enable airflow in DEM config
	sim.config.enable_airflow = true;

	// Initialize DEM
	sim.init();

	// Compute DEM substeps per LBM step
	// DEM dt is set by sim.init() (auto-dt), we need enough substeps to cover dt_lbm_si
	if(sim.config.dt > 0.0 && dt_lbm_si > 0.0f) {
		int substeps = (int)std::ceil((double)dt_lbm_si / sim.config.dt);
		if(substeps > 1) {
			sim.config.dt = (double)dt_lbm_si / (double)substeps;
		}
		std::cout << "[DEM-LBM] DEM substeps per LBM step: " << substeps
		          << " (DEM dt=" << sim.config.dt << " s, LBM dt=" << dt_lbm_si << " s)\n";
	}

	std::cout << "[DEM-LBM] Coupling initialized: grid=" << Nx << "x" << Ny << "x" << Nz
	          << " spacing=" << spacing*1000.0f << " mm"
	          << " u_conv=" << u_conversion << " m/s per lu"
	          << " f_conv=" << f_conversion
	          << (twoway ? " TWO-WAY" : " ONE-WAY") << "\n";
}

void DemCoupling::step(LBM& lbm) {
	const ulong N = (ulong)Nx * (ulong)Ny * (ulong)Nz;

	// 1. Read LBM velocity field from GPU to CPU
	lbm.u.read_from_device();

	// 2. Convert LBM velocity to SI and store in buffers
	for(ulong n = 0ull; n < N; n++) {
		ux_buf[n] = lbm.u.x[n] * u_conversion;
		uy_buf[n] = lbm.u.y[n] * u_conversion;
		uz_buf[n] = lbm.u.z[n] * u_conversion;
	}

	// 3. Update DEM airflow field with new velocity data
	sim.airflow.update_from_arrays(ux_buf.data(), uy_buf.data(), uz_buf.data());

	// 4. Step DEM forward by one LBM timestep (multiple DEM substeps)
	//    Accumulate drag forces from all substeps for two-way feedback
	std::vector<dem::AirflowField::ParticleDrag> accumulated_drags;
	int substeps = std::max(1, (int)std::ceil((double)dt_lbm_si / sim.config.dt));
	for(int i = 0; i < substeps; i++) {
		if(sim.is_finished()) break;
		sim.step();
		if(twoway) {
			// Collect drag forces from this substep (already computed in sim.compute_forces())
			for(const auto& d : sim.last_drag_forces) {
				accumulated_drags.push_back(d);
			}
		}
	}

	// 5. Two-way coupling: distribute reaction forces back to LBM force field
#ifdef FORCE_FIELD
	if(twoway && !accumulated_drags.empty()) {
		// Distribute DEM drag reaction forces onto grid (SI force-per-volume)
		sim.airflow.distribute_forces_to_grid(accumulated_drags, fx_buf.data(), fy_buf.data(), fz_buf.data());

		// Scale by 1/substeps to get time-averaged force per LBM step,
		// then convert SI force-per-volume to LBM lattice force units
		const float scale = f_conversion / (float)substeps;

		// Write force field into LBM F arrays (CPU side)
		for(ulong n = 0ull; n < N; n++) {
			lbm.F.x[n] = fx_buf[n] * scale;
			lbm.F.y[n] = fy_buf[n] * scale;
			lbm.F.z[n] = fz_buf[n] * scale;
		}

		// Upload force field to GPU
		lbm.F.write_to_device();
	}
#endif // FORCE_FIELD

	// 6. Update visualization buffer
	update_vis();
}

void DemCoupling::update_vis() {
	std::lock_guard<std::mutex> lock(vis_mutex);
	vis_particles.clear();

	const float inv_spacing = 1.0f / spacing;
	for(const auto& p : sim.particles.particles) {
		if(!p.active) continue;
		const auto& st = sim.particles.type_of(p);

		DemParticleVis pv;
		// Convert DEM position (SI meters) to LBM centered coordinates (lattice units)
		pv.x = (float)(p.pos.x * (double)inv_spacing);
		pv.y = (float)(p.pos.y * (double)inv_spacing);
		pv.z = (float)(p.pos.z * (double)inv_spacing);

		// Convert radius to lattice units
		double r = (st.shape == dem::SeedShape::SPHERE) ? st.radius : std::max({st.a, st.b, st.c});
		pv.radius = (float)(r * (double)inv_spacing);

		// Color based on velocity magnitude
		double v = p.vel.length();
		float vnorm = (float)std::min(v / 2.0, 1.0); // normalize to [0,1] for 0..2 m/s
		int ri = (int)(255.0f * vnorm);
		int bi = (int)(255.0f * (1.0f - vnorm));
		int gi = (int)(128.0f * (1.0f - std::fabs(vnorm - 0.5f) * 2.0f));
		pv.color = (ri << 16) | (gi << 8) | bi;

		vis_particles.push_back(pv);
	}
}

#ifdef GRAPHICS
void DemCoupling::draw_particles() {
	std::lock_guard<std::mutex> lock(vis_mutex);
	for(const auto& pv : vis_particles) {
		float3 pos(pv.x, pv.y, pv.z);
		// Draw filled circle: draw concentric circles from large to small
		float r = std::max(pv.radius, 1.0f); // minimum 1 lattice unit for visibility
		draw_circle(pos, r, pv.color);
		if(r > 1.5f) draw_circle(pos, r * 0.7f, pv.color);
		if(r > 2.5f) draw_circle(pos, r * 0.4f, pv.color);
		draw_pixel(pos, pv.color); // always draw center pixel
	}
}
#else
void DemCoupling::draw_particles() {} // no-op without graphics
#endif // GRAPHICS

#endif // DEM_COUPLING
