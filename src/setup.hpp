#pragma once

#include "defines.hpp"
#include "lbm.hpp"
#include "shapes.hpp"

#ifdef DEM_COUPLING
#include "dem_coupling.hpp"
#endif // DEM_COUPLING

void main_setup(); // main setup script