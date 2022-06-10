//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>

namespace py = pybind11;

extern "C" {

void IncCVMSWells_copy_states_from_reservoir();
void MSWellsData_init();
void IncPrimSecdMSWells_compute();
void LoisThermoHydroMSWells_compute();
void LeafMSWells_init_data();
void VSHydroMSWells_init();
void VSHydroMSWells_compute();
}

struct MSWell {};

void add_mswell_wrappers(py::module& module) {
   py::class_<MSWell>(module, "MSWell").def(py::init<>());

   module.def("mswells_copy_states_from_reservoir",
              &IncCVMSWells_copy_states_from_reservoir,
              "For each vertex of all  multi-segmented producer wells, copy "
              "the  Coats variables from the reservoir");

   module.def("mswells_init_edge_data", &MSWellsData_init);
   module.def("IncPrimSecdMSWells_compute", &IncPrimSecdMSWells_compute);
   module.def("LoisThermoHydroMSWells_compute",
              &LoisThermoHydroMSWells_compute);

   module.def("mswells_init_leaf_data", &LeafMSWells_init_data);

   module.def("VSHydroMSWells_init", &VSHydroMSWells_init);
   module.def("VSHydroMSWells_compute", &VSHydroMSWells_compute);
}
