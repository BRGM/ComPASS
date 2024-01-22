//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <metis.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <iostream>
#include <limits>
// #include <span>
// #include <type_traits>

// #include "Metis_wrapper.h"

namespace py = pybind11;

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
   //    std::cerr << "METIS: nb_vertices = " << nb_vertices << std::endl;
   //    std::cerr << "METIS: nb_parts = " << nb_parts << std::endl;
   //    auto offprox = offsets.template mutable_unchecked<1>();
   //    std::cerr << "METIS: offsets = ";
   //    for (auto&& i : std::span(offsets.data(), (std::size_t)(nb_vertices +
   //    1))) {
   //       std::cerr << " " << i;
   //    }
   //    std::cerr << std::endl;
   //    std::cerr << "METIS: neighbors = ";
   //    for (auto&& i :
   //         std::span(neighbors.data(), (std::size_t)(offprox(nb_vertices))))
   //         {
   //       std::cerr << " " << i;
   //    }
   //    std::cerr << std::endl;
   //    idx_t* res = new idx_t[nb_vertices];
   auto status =
       method(&nb_vertices, &nb_constraints, offsets.mutable_data(),
              neighbors.mutable_data(), nullptr, nullptr, nullptr, &nb_parts,
              nullptr, nullptr, options, &edge_cut, result.mutable_data());
   //    auto status = method(&nb_vertices, &nb_constraints,
   //    offsets.mutable_data(),
   //                         neighbors.mutable_data(), nullptr, nullptr,
   //                         nullptr, &nb_parts, nullptr, nullptr, options,
   //                         &edge_cut, res);

   assert(status == METIS_OK);

   //    std::cerr << "METIS: result = ";
   //    for (auto&& i : std::span(res, (std::size_t)nb_vertices)) {
   //       std::cerr << " " << i;
   //    }
   //    std::cerr << std::endl;

   //    std::ranges::copy(std::span(res, (std::size_t)nb_vertices),
   //                      result.mutable_data());

   //    delete[] res;

   return result;
}

PYBIND11_MODULE(metis, module) {
   module.doc() =
       "quick and dirty wrapper to metis partitioning used in ComPASS";
   module.def("part_graph", &part_graph,
              "call Metis part graph methods, defaults to recursive algorithm",
              py::arg("neighbors"), py::arg("offsets"), py::arg("nb_parts"),
              py::arg("Kway") = false);
}
