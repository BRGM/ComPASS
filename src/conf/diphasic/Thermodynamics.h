#pragma once

#include "../common/fortran_thermodynamics.h"

extern "C" {
void DiphasicFlash_enforce_consistent_molar_fractions(X &);
}
