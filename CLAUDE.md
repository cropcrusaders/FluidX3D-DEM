# CLAUDE.md

Guide for AI assistants working on this codebase.

## Project Overview

**FluidX3D-DEM** combines two simulation components:

- **FluidX3D** (`src/`): GPU-accelerated lattice Boltzmann method (LBM) CFD solver using OpenCL. Created by Dr. Moritz Lehmann. Runs on all GPUs (NVIDIA, AMD, Intel) and CPUs via OpenCL.
- **DEM module** (`dem/`): Discrete Element Method particle simulator for seed drop-tube analysis. Simulates rigid-body particle dynamics with collision detection, contact forces, and optional one-way fluid coupling from FluidX3D velocity fields.

The two components are loosely coupled: FluidX3D exports a 3D velocity field, which the DEM module loads and interpolates to compute drag forces on particles.

## Repository Structure

```
FluidX3D-DEM/
├── src/                          # FluidX3D core (C++17, OpenCL)
│   ├── main.cpp                  # Entry point
│   ├── setup.cpp                 # Simulation setups (uncomment to select)
│   ├── lbm.hpp/cpp               # LBM_Domain & LBM classes, GPU memory, multi-GPU
│   ├── kernel.hpp/cpp            # OpenCL kernel compilation & management
│   ├── graphics.hpp/cpp          # Rendering (rasterization, raytracing)
│   ├── opencl.hpp                # OpenCL C++ wrapper
│   ├── utilities.hpp             # Vector math, threading, file I/O, memory pools
│   ├── units.hpp                 # SI unit conversion (LBM <-> physical)
│   ├── shapes.hpp/cpp            # Geometric primitives, STL voxelization
│   ├── info.hpp/cpp              # Runtime status reporting
│   ├── defines.hpp               # Compile-time feature flags (velocity set, extensions)
│   ├── lodepng.hpp/cpp           # Embedded PNG library
│   ├── OpenCL/                   # Vendored OpenCL headers and libs
│   └── X11/                      # Vendored X11 headers and libs (Linux interactive)
├── dem/                          # DEM seed drop-tube simulator (C++17, CMake)
│   ├── include/                  # Headers
│   │   ├── simulator.hpp         # Main simulation orchestrator
│   │   ├── particle.hpp          # Particle state, SeedType, material properties
│   │   ├── collision.hpp         # Spatial hash grid broadphase + narrowphase
│   │   ├── contact.hpp           # Spring-dashpot contact model with friction
│   │   ├── geometry.hpp          # Triangle mesh (BVH) + analytic primitives
│   │   ├── airflow.hpp           # Velocity field interpolation, drag computation
│   │   ├── injector.hpp          # Particle source configuration
│   │   ├── metrics.hpp           # Exit detection, per-seed statistics
│   │   ├── dem_math.hpp          # vec3, quaternion, mat3, AABB, RNG
│   │   └── config_parser.hpp     # JSON config parser (no external deps)
│   ├── src/                      # Implementations (one .cpp per header)
│   ├── tests/test_main.cpp       # 14 unit tests, 53 checks
│   ├── examples/                 # JSON config templates (6 configs)
│   ├── python/dem_postprocess.py # Matplotlib postprocessing
│   └── CMakeLists.txt            # CMake build (dem_core lib, dem_sim, dem_tests)
├── makefile                      # GNU Make for FluidX3D
├── make.sh                       # Platform-detection build script for FluidX3D
├── FluidX3D.sln / .vcxproj      # Visual Studio project (Windows)
├── README.md                     # FluidX3D documentation
├── DOCUMENTATION.md              # Getting started guide
└── LICENSE.md                    # Custom license (non-commercial, see below)
```

## Build Commands

### FluidX3D (GNU Make + OpenCL)

```bash
# Auto-detect platform and build (also runs the binary on success)
./make.sh

# Build only (without running)
make Linux-X11 -j$(nproc)   # Linux with X11 interactive graphics
make Linux -j$(nproc)        # Linux without X11
make macOS -j$(nproc)        # macOS
make Android -j$(nproc)      # Android (Termux)

# Clean
make clean
```

- Compiler: `g++` with `-std=c++17 -pthread -O -Wno-comment`
- Output: `bin/FluidX3D`
- Intermediate objects: `temp/`
- Requires OpenCL runtime (GPU drivers or CPU fallback)
- X11 needed only for `INTERACTIVE_GRAPHICS` mode on Linux

### DEM Module (CMake)

```bash
cd dem
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

- Compiler flags: `-O3 -Wall -Wextra -Wpedantic` (Release)
- Produces: `dem_sim` (CLI), `dem_tests` (unit tests)
- Optional dependency: OpenMP (auto-detected)
- No external C++ libraries required

### Running DEM Tests

```bash
cd dem/build
./dem_tests          # runs all 14 tests (53 checks)
# or via CTest:
ctest
```

Tests cover: vector math, quaternion rotation, AABB operations, triangle closest-point, RNG determinism, JSON parsing, seed inertia, drag force, Schiller-Naumann drag coefficient, cylinder collision, free-fall trajectory, bounce restitution, inclined-plane friction, and elastic energy conservation.

### Running DEM Simulations

```bash
cd dem/build
./dem_sim ../examples/straight_tube.json
./dem_sim ../examples/curved_tube.json
./dem_sim ../examples/spiral_tube.json

# With CLI overrides:
./dem_sim config.json --output-dir output/ --dt 1e-5 --max-time 2.0 --seed 12345
```

## Configuration

### FluidX3D (`src/defines.hpp`)

Configuration is compile-time via preprocessor macros. Key settings:

- **Velocity set**: Uncomment one of `D2Q9`, `D3Q15`, `D3Q19` (default), `D3Q27`
- **Collision operator**: `SRT` (default) or `TRT`
- **Precision**: `FP16S` (default, range-shifted FP16) or `FP16C` or neither (FP32)
- **Extensions**: `VOLUME_FORCE`, `FORCE_FIELD`, `EQUILIBRIUM_BOUNDARIES`, `MOVING_BOUNDARIES`, `SURFACE`, `TEMPERATURE`, `SUBGRID`, `PARTICLES`
- **Graphics**: `INTERACTIVE_GRAPHICS`, `INTERACTIVE_GRAPHICS_ASCII`, `GRAPHICS`
- **Benchmark**: `BENCHMARK` (disables all extensions, runs benchmark setup)

Simulation setups are in `src/setup.cpp` - uncomment the desired `main_setup()` block.

### DEM Module (JSON configs)

Configuration is runtime via JSON files. See `dem/examples/` for templates. Key sections:

- `simulation`: timestep, max_time, gravity, RNG seed, output settings
- `seed_types`: particle shapes (sphere/ellipsoid), mass, radius, material
- `contact`: spring stiffness (kn, kt), restitution, friction coefficients
- `geometry`: triangle meshes (STL/OBJ) and analytic primitives (cylinder, cone, plane)
- `injector`: particle source position, rate, aperture shape
- `airflow`: enable/disable fluid coupling, velocity field file path, drag model

## Architecture

### FluidX3D Key Classes

- **`LBM`**: Multi-GPU orchestrator. Manages domain decomposition and global state.
- **`LBM_Domain`**: Single-GPU domain. Owns GPU memory, kernels, and halo exchange.
- **`Device`**: OpenCL device wrapper (context, queues, memory info).
- **`Kernel`**: OpenCL kernel with compile-time options injection.
- **`Memory<T>`**: GPU/CPU buffer pair with async transfer (RAII).
- **`Graphics`**: Frame rendering pipeline (rasterization and raytracing via OpenCL kernels).

### DEM Key Classes (all in `namespace dem`)

- **`Simulator`**: Time-stepping orchestrator. Owns all subsystems.
- **`Particle`**: 6-DOF state (position, velocity, quaternion, angular velocity, forces/torques).
- **`CollisionDetector`**: Spatial hash grid broadphase + sphere narrowphase.
- **`ContactModel`**: Linear spring-dashpot with Coulomb friction and rolling friction.
- **`Geometry`**: Collection of triangle meshes (BVH-accelerated) and analytic primitives.
- **`AirflowField`**: Loads FluidX3D velocity field, trilinear interpolation, Schiller-Naumann drag.
- **`Injector`**: Configurable particle source.
- **`MetricsCollector`**: Exit detection, per-seed statistics, CSV output.

### Data Flow (Fluid-Particle Coupling)

```
FluidX3D (GPU, OpenCL)              DEM (CPU)
    │                                   │
    ├─ LBM timesteps ──► Export         │
    │   velocity field (binary/NPY)     │
    │                                   │
    │                     Load ◄────────┤
    │                     Trilinear     │
    │                     interpolation │
    │                         │         │
    │                     Drag force ───┤
    │                     (Schiller-    │
    │                      Naumann)     │
    │                                   │
    │                     Output ◄──────┤
    │                     (CSV, VTK)    │
```

Currently one-way coupling only (DEM reads FluidX3D output; does not feed back).

## Coding Conventions

### General

- **C++ standard**: C++17 throughout
- **No heavy external dependencies**: Math, JSON parsing, PNG I/O are all self-contained
- **One concept per file** in the DEM module (header + implementation)

### FluidX3D Style

- Classes: `CamelCase` with underscores for compound names (`LBM_Domain`)
- Functions/variables: `snake_case`
- Constants/macros: `UPPERCASE` (`TYPE_S`, `VIS_FLAG_LATTICE`)
- Heavy use of preprocessor macros for compile-time configuration
- Single large headers (e.g., `utilities.hpp` contains vector math, threading, file I/O)
- Minimal comments; self-documenting via clear naming
- Memory layout: Structure-of-arrays (SoA) for GPU cache efficiency
- Linear indexing: `x + Nx*(y + Ny*z)`

### DEM Style

- All code in `namespace dem`
- Classes/structs: `CamelCase` (`Particle`, `SeedType`, `ContactProperties`)
- Variables: `snake_case` (`pos`, `vel`, `omega`, `dt`)
- No macros (except include guards and `#ifdef _OPENMP`)
- Extensive comments explaining physics decisions
- Custom test framework using `CHECK()` / `CHECK_NEAR()` macros

### Things to Avoid

- Do not add external library dependencies (Eigen, Boost, nlohmann/json, etc.)
- Do not change `defines.hpp` flags without understanding the compile-time cascade (e.g., `SURFACE` forces `UPDATE_FIELDS`)
- Do not modify vendored code in `src/OpenCL/`, `src/X11/`, or `src/lodepng.*`
- Be aware of the `BENCHMARK` macro which disables all extensions and graphics

## License Constraints

Custom proprietary license (see `LICENSE.md`). Key restrictions:
- No commercial use
- No military use
- No AI model training on source code
- Altered versions must be published if binaries/results are published
- Must cite referenced papers in scientific publications

## Key File Locations

| Purpose | Path |
|---------|------|
| FluidX3D entry point | `src/main.cpp` |
| Simulation setup selection | `src/setup.cpp` |
| Compile-time feature flags | `src/defines.hpp` |
| FluidX3D build script | `make.sh` |
| DEM entry point | `dem/src/main.cpp` |
| DEM build config | `dem/CMakeLists.txt` |
| DEM unit tests | `dem/tests/test_main.cpp` |
| DEM example configs | `dem/examples/*.json` |
| DEM postprocessing | `dem/python/dem_postprocess.py` |
