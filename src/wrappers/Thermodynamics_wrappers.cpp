//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

#include <cstdint>

#include "StateObjects.h"

namespace py = pybind11;

using Physical_property_with_derivatives_ptr = void (*)(const Xalpha &,
                                                        double *, Xalpha &);
using Physical_property_ptr = double (*)(const Xalpha &);
Physical_property_with_derivatives_ptr phy_prop_d[ComPASS_NUMBER_OF_PHASES];
Physical_property_ptr phy_prop[ComPASS_NUMBER_OF_PHASES];

extern "C" {
void register_c_viscosities_with_derivatives(decltype(phy_prop_d));
void register_c_viscosities_without_derivatives(decltype(phy_prop));
}

void add_physical_properties_wrappers(py::module &module) {
   module.def("register_c_viscosities_with_derivatives", [](py::args args) {
      if (py::len(args) != ComPASS_NUMBER_OF_PHASES)
         throw std::runtime_error("wrong number of phases");
      for (std::size_t k = 0; k < ComPASS_NUMBER_OF_PHASES; ++k) {
         phy_prop_d[k] =
             reinterpret_cast<Physical_property_with_derivatives_ptr>(
                 args[k].cast<std::uintptr_t>());
      }
      register_c_viscosities_with_derivatives(phy_prop_d);
   });
   module.def("register_c_viscosities_without_derivatives", [](py::args args) {
      if (py::len(args) != ComPASS_NUMBER_OF_PHASES)
         throw std::runtime_error("wrong number of phases");
      for (std::size_t k = 0; k < ComPASS_NUMBER_OF_PHASES; ++k) {
         phy_prop[k] = reinterpret_cast<Physical_property_ptr>(
             args[k].cast<std::uintptr_t>());
      }
      register_c_viscosities_without_derivatives(phy_prop);
   });
}
