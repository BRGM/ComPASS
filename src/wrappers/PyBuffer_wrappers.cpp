//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "PyBuffer_wrappers.h"

void add_pybuffer_wrappers(py::module& module) {
   IntBuffer::add_buffer_class(module, "IntBuffer");
   DoubleBuffer::add_buffer_class(module, "DoubleBuffer");
   CoordinatesBuffer::add_buffer_class(module, "CoordinatesBuffer");
   TensorBuffer::add_buffer_class(module, "TensorBuffer");
}
