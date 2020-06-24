//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "Metis_wrapper.h"

#include <metis.h>
#include <pybind11/numpy.h>

#include <cassert>
#include <type_traits>

#include "COC.h"

auto part_graph(COC& adjacencies, idx_t nb_parts, const bool Kway) {
   static_assert(std::is_same<idx_t, COC::value_type>::value,
                 "unconsistent types");

   auto nb_vertices = static_cast<idx_t>(adjacencies.number_of_containers());
   auto nb_constraints = static_cast<idx_t>(1);
   idx_t edge_cut;
   auto result = py::array_t<idx_t, py::array::c_style>{
       adjacencies.number_of_containers()};

   idx_t options[METIS_NOPTIONS];
   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_NUMBERING] = 1;  // Fortran-style numbering

   auto method = METIS_PartGraphRecursive;
   if (Kway) method = METIS_PartGraphKway;
   auto status = method(
       &nb_vertices, &nb_constraints, adjacencies.mutable_offset_data(),
       adjacencies.mutable_content_data(), nullptr, nullptr, nullptr, &nb_parts,
       nullptr, nullptr, options, &edge_cut, result.mutable_data());

   assert(status == METIS_OK);

   return result;
}

void add_Metis_wrapper(py::module& module) {
   module.def("metis_part_graph", &part_graph,
              "call Metis part graph methods, defaults to recursive algorithm",
              py::arg("adjacencies"), py::arg("nb_parts"),
              py::arg("Kway") = false);
}
