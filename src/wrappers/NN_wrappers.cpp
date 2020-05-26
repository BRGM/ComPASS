//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "ArrayWrapper.h"
#include "StringWrapper.h"

// Fortran functions
extern "C" {
void NN_init_warmup(const StringWrapper&);
void NN_wait_for_debug(bool);
void NN_main_summarize_timestep();
void NN_finalize();
void NN_init_phase2_partition(ArrayWrapper&);
void init_phase2_build_local_mesh();
void init_phase2_setup_contexts();
void init_phase2_setup_solvers();
}

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

#include "NN_wrappers.h"

void add_NN_wrappers(py::module& module) {
   module.def(
       "init_warmup",
       [](const std::string& LogFile) { NN_init_warmup(LogFile); },
       "Initialisation of ComPASS - warmup phase.");

   module.def(
       "init_phase2_partition",
       [](py::array_t<int, py::array::c_style> colors) {
          assert(colors.ndim() == 1);
          auto wrapper =
              ArrayWrapper::wrap(colors.mutable_data(), colors.size());
          NN_init_phase2_partition(wrapper);
       },
       "Partition mesh.");

   module.def("init_phase2_build_local_mesh", &init_phase2_build_local_mesh);
   module.def("init_phase2_setup_contexts", &init_phase2_setup_contexts);
   module.def("init_phase2_setup_solvers", &init_phase2_setup_solvers);

   module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");
   module.def("init_wait_for_debug", &NN_wait_for_debug);
}
