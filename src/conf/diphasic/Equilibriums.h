#pragma once

#include "StateObjects.h"

using Real = X::Model::Real;
using Phase_vector = X::Model::Phase_vector;

Phase_vector diphasic_equilibrium(const Phase_vector &pa, const Real &T,
                                  const double atol = 1e-7,
                                  std::size_t maxiter = 1000);
