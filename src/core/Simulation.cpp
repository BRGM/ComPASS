//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "Simulation.h"

namespace ComPASS {

/** Alternative 1: the simulation object is on the C++ side *

std::unique_ptr<Simulation> simulation_handle;
py::object simulation_pyhandle = py::none();

void create_simulation_object(int proc) {
   if (simulation_handle) {
      assert(!simulation_pyhandle.is_none());
      throw std::runtime_error("A simulation object has already been defined.");
   }
   assert(simulation_pyhandle.is_none());
   simulation_handle = std::make_unique<Simulation>(proc);
   simulation_pyhandle = py::cast(simulation_handle.get());
   assert(simulation_handle.get() == simulation_pyhandle.cast<Simulation*>());
}

Simulation* get_simulation_handle() { return simulation_handle.get(); }
py::object& get_simulation_pyhandle() { return simulation_pyhandle; }

/**/

/** Alternative 2: Use the python simulation object */
py::object simulation = py::none();

void register_simulation(py::object& obj) {
   if (!simulation.is_none())
      throw std::runtime_error(
          "A simulation object has already been registered.");
   simulation = obj;
   assert(!simulation.is_none());
}

void release_simulation() {
   if (!simulation.is_none()) simulation = py::none();
}

py::object& get_simulation() { return simulation; }

/**/

}  // namespace ComPASS

void add_simulation_pyinternals(py::module& module) {
   /** Alternative 1: the simulation object is on the C++ side *

   py::class_<ComPASS::Simulation>(module, "Simulation");
   module.def("create_simulation_object", &ComPASS::create_simulation_object,
              py::arg("proc") = -1);

   /**/

   /** Alternative 2: Use the python simulation object */
   module.def("register_simulation", &ComPASS::register_simulation);
   module.def("release_simulation", &ComPASS::release_simulation);
}
