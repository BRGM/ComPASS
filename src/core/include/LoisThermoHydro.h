#pragma once

#include "StateObjects.h"

using Phase_vector = X::Model::Phase_vector;

extern "C" {

void fill_kr_arrays(const std::size_t, const int, X*, int*, double*, double*);

void fill_phase_pressure_arrays(const std::size_t, const int, X*, int*, double*,
                                double*);
}

Phase_vector phase_pressures(X&);
