//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "Simulation.h"

#include <memory>

namespace ComPASS {

std::unique_ptr<Simulation> simulation_handle;

void create_simulation_object(int proc) {
   if (simulation_handle)
      throw std::runtime_error("A simulation object has already been defined.");
   simulation_handle = std::make_unique<Simulation>(proc);
}

Simulation* simulation_object() { return simulation_handle.get(); }

}  // namespace ComPASS

#include <pybind11/pybind11.h>
namespace py = pybind11;

void add_simulation_pyinternals(py::module& module) {
   module.def("create_simulation_object", &ComPASS::create_simulation_object,
              py::arg("proc") = -1);
}
