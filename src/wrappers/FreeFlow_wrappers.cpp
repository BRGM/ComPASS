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

#include "ArrayWrapper.h"

#ifdef _WIP_FREEFLOW_STRUCTURES_
// Fortran functions
extern "C" {
void clear_freeflow_faces();
void set_freeflow_faces(const ArrayWrapper &);
}
#endif
namespace py = pybind11;

void add_freeflow_wrappers(py::module &module) {
#ifdef _WIP_FREEFLOW_STRUCTURES_
   module.attr("has_freeflow_structures") = true;
#else
   module.attr("has_freeflow_structures") = false;
#endif

#ifdef _WIP_FREEFLOW_STRUCTURES_
   module.def("clear_freeflow_faces", &clear_freeflow_faces);
   module.def(
       "set_freeflow_faces",
       [](py::array_t<int, py::array::c_style | py::array::forcecast> faces) {
          if (faces.ndim() != 1)
             throw std::runtime_error(
                 "Faces array shoud be a single dimension array.");
          auto wrapper =
              ArrayWrapper::wrap(faces.mutable_data(0), faces.size());
          set_freeflow_faces(wrapper);
       });
#endif
}
