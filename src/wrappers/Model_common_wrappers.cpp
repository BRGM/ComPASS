//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>

#include "Model_wrappers.h"

// Fortran functions
extern "C" {
int model_number_of_phases();
int model_number_of_components();
int model_number_of_contexts();
bool is_context_locked(int);
void lock_context(int);
void unlock_context(int);
}

void add_common_model_wrappers(py::module& module) {
   module.def("init_model", &init_model);
   module.def("finalize_model", &finalize_model);
   module.def("model_number_of_phases", &model_number_of_phases);
   module.def("model_number_of_components", &model_number_of_components);
   module.def("model_number_of_contexts", &model_number_of_contexts);
}

void use_context_locker(py::module& module) {
   for (int k = 1; k <= model_number_of_contexts(); ++k) {
      unlock_context(k);
      assert(!is_context_locked(k));
   }
   module.def("is_context_locked", &is_context_locked);
   module.def("lock_context", &lock_context);
   module.def("unlock_context", &unlock_context);
}

void add_model_wrappers(py::module& module) {
   add_common_model_wrappers(module);
   add_specific_model_wrappers(module);
}
