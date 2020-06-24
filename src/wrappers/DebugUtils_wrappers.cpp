//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

// Fortran functions
extern "C" {
void debug_utils_dump_mesh_info();
}

#include "DebugUtils_wrappers.h"

void add_debug_utils_wrappers(py::module& module) {
   module.def("debug_utils_dump_mesh_info", &debug_utils_dump_mesh_info);
}
