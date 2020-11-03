#include <metis.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <type_traits>
namespace py = pybind11;

auto part_graph(
    py::array_t<idx_t, py::array::c_style | py::array::forcecast>& neighbors,
    py::array_t<idx_t, py::array::c_style | py::array::forcecast>& offsets,
    idx_t nb_parts, const bool Kway) {
   auto nb_vertices = static_cast<idx_t>(offsets.size() - 1);
   auto nb_constraints = static_cast<idx_t>(1);
   idx_t edge_cut = -1;
   auto result = py::array_t<idx_t, py::array::c_style>{
       static_cast<std::size_t>(nb_vertices)};

   idx_t options[METIS_NOPTIONS];
   METIS_SetDefaultOptions(options);
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

PYBIND11_MODULE(metis, module) {
   module.doc() =
       "quick and dirty wrapper to metis partitioning used in ComPASS";
   module.def("metis_part_graph", &part_graph,
              "call Metis part graph methods, defaults to recursive algorithm",
              py::arg("neighbors"), py::arg("offsets"), py::arg("nb_parts"),
              py::arg("Kway") = false);
}
