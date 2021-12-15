//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

#include <cassert>

#include "NewtonIncrements.h"

// Fortran functions
extern "C" {
void Jacobian_JacBigA_BigSm(const double&);
void Jacobian_Regularization();
void Jacobian_Schur();
void Jacobian_Alignment_man();
void Jacobian_Alignment_diag();
void Jacobian_GetSolCell(NewtonIncrements::Pointers<double>);
}

namespace py = pybind11;

void add_Jacobian_wrappers(py::module& module) {
   module.def("Jacobian_JacBigA_BigSm", &Jacobian_JacBigA_BigSm);
   module.def("Jacobian_Regularization", &Jacobian_Regularization);
   module.def("Jacobian_Schur", &Jacobian_Schur);
   module.def("Jacobian_Alignment_man", &Jacobian_Alignment_man);
   module.def("Jacobian_Alignment_diag", &Jacobian_Alignment_diag);
   module.def("Jacobian_GetSolCell", [](NewtonIncrements& increments) {
      //   py::print("sizes:", nb_primary_variables(), nb_nodes(),
      //       		    nb_fractures(), nb_cells(), nb_injectors(),
      //       nb_producers());
      Jacobian_GetSolCell(increments.pointers());
   });
}
