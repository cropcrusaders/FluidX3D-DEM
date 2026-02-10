cmake_minimum_required(VERSION 3.14)
project(dem_seed_sim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

# OpenMP (optional but recommended)
find_package(OpenMP)

# Source files
set(DEM_SOURCES
    src/geometry.cpp
    src/contact.cpp
    src/collision.cpp
    src/airflow.cpp
    src/injector.cpp
    src/metrics.cpp
    src/simulator.cpp
    src/config_parser.cpp
)

# Static library
add_library(dem_core STATIC ${DEM_SOURCES})
target_include_directories(dem_core PUBLIC include)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dem_core PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(dem_core PRIVATE -Wall -Wextra -Wpedantic)

# CLI executable
add_executable(dem_sim src/main.cpp)
target_link_libraries(dem_sim PRIVATE dem_core)

# Tests
enable_testing()
add_executable(dem_tests tests/test_main.cpp)
target_link_libraries(dem_tests PRIVATE dem_core)
target_include_directories(dem_tests PRIVATE include)
add_test(NAME dem_unit_tests COMMAND dem_tests)

# Install
install(TARGETS dem_sim DESTINATION bin)
install(DIRECTORY examples/ DESTINATION share/dem_seed_sim/examples)
