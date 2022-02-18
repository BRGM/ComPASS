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
#include <limits>
#include <type_traits>

#include "COC.h"

auto part_graph(
    py::array_t<idx_t, py::array::c_style | py::array::forcecast>& neighbors,
    py::array_t<idx_t, py::array::c_style | py::array::forcecast>& offsets,
    idx_t nb_parts, const bool Kway) {
   assert(offsets.size() > 0);
   auto nb_vertices = static_cast<idx_t>(offsets.size() - 1);
#ifndef NDEBUG
   assert(neighbors.ndim() == 1);
   auto pk = neighbors.data(0);
   for (std::size_t k = 0; k < neighbors.size(); ++k) {
      if (*pk < 0 || *pk >= nb_vertices) {
         throw std::runtime_error(
             "Indexing of neighbors should follow C-style.");
      }
      ++pk;
   }
#endif
   auto nb_constraints = static_cast<idx_t>(1);
   idx_t edge_cut;
   static_assert(std::numeric_limits<decltype(nb_vertices)>::max() <
                 std::numeric_limits<py::ssize_t>::max());
   auto result = py::array_t<idx_t, py::array::c_style>{
       static_cast<py::ssize_t>(nb_vertices)};

   idx_t options[METIS_NOPTIONS];
   METIS_SetDefaultOptions(options);
   // We will use the default seed for Kway part graph
   // but you can set yout own with: options[METIS_OPTION_SEED] = *your seed
   // here* WARNING: if you set the Fortran-style numbering
   //          It seems that offsets must also be start at one
   //          (which seems weird...)
   // options[METIS_OPTION_NUMBERING] = 1;  // Fortran-style numbering

   auto method = METIS_PartGraphRecursive;
   if (Kway) method = METIS_PartGraphKway;
   auto status =
       method(&nb_vertices, &nb_constraints, offsets.mutable_data(),
              neighbors.mutable_data(), nullptr, nullptr, nullptr, &nb_parts,
              nullptr, nullptr, options, &edge_cut, result.mutable_data());

   assert(status == METIS_OK);

   return result;
}

void add_Metis_wrapper(py::module& module) {
   module.def("metis_part_graph", &part_graph,
              "call Metis part graph methods, defaults to recursive algorithm",
              py::arg("neighbors"), py::arg("offsets"), py::arg("nb_parts"),
              py::arg("Kway") = false);
}
