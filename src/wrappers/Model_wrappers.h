//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_model();
void finalize_model();
void add_model_wrappers(py::module& module);
void add_common_model_wrappers(py::module& module);
void add_specific_model_wrappers(py::module& module);
void use_context_locker(py::module& module);
