#pragma once

#include "../common/fortran_thermodynamics.h"

extern "C" {
void FluidThermodynamics_coeff_CO2_liquid_fug(const double &, const double &,
                                              double &, double &, double &);
}
