#pragma once

#include "StateObjects.h"

extern "C" {
void fill_kr_arrays(const std::size_t, const int, X*, int*, double*, double*);
void fill_phase_pressure_arrays(const std::size_t, const int, X*, int*,
                                double*);
}

void update_phase_pressures(X&);
