cmake_minimum_required(VERSION 3.14)
project(coupled_sim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

# OpenMP (optional)
find_package(OpenMP)

# ─── DEM library (from sibling directory) ──────────────────────────────────
set(DEM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../dem)
set(DEM_SOURCES
    ${DEM_DIR}/src/geometry.cpp
    ${DEM_DIR}/src/contact.cpp
    ${DEM_DIR}/src/collision.cpp
    ${DEM_DIR}/src/airflow.cpp
    ${DEM_DIR}/src/injector.cpp
    ${DEM_DIR}/src/metrics.cpp
    ${DEM_DIR}/src/simulator.cpp
    ${DEM_DIR}/src/config_parser.cpp
)

add_library(dem_core STATIC ${DEM_SOURCES})
target_include_directories(dem_core PUBLIC ${DEM_DIR}/include)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dem_core PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(dem_core PRIVATE -Wall -Wextra -Wpedantic)

# ─── FluidX3D sources (compiled as library, requires OpenCL) ──────────────
set(LBM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(LBM_SOURCES
    ${LBM_DIR}/lbm.cpp
    ${LBM_DIR}/kernel.cpp
    ${LBM_DIR}/graphics.cpp
    ${LBM_DIR}/info.cpp
    ${LBM_DIR}/shapes.cpp
    ${LBM_DIR}/lodepng.cpp
)

add_library(lbm_core STATIC ${LBM_SOURCES})
target_include_directories(lbm_core PUBLIC ${LBM_DIR} ${LBM_DIR}/OpenCL/include)
target_compile_options(lbm_core PRIVATE -std=c++17 -pthread -O -Wno-comment)

# OpenCL linking (platform-dependent)
if(APPLE)
    target_link_libraries(lbm_core PUBLIC "-framework OpenCL")
else()
    target_link_directories(lbm_core PUBLIC ${LBM_DIR}/OpenCL/lib)
    target_link_libraries(lbm_core PUBLIC OpenCL pthread)
endif()

# ─── Coupled executable ───────────────────────────────────────────────────
add_executable(coupled_sim coupled_sim.cpp)
target_link_libraries(coupled_sim PRIVATE dem_core lbm_core)
target_include_directories(coupled_sim PRIVATE ${LBM_DIR} ${DEM_DIR}/include)

install(TARGETS coupled_sim DESTINATION bin)
