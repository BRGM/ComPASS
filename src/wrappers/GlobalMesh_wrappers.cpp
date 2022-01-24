//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>

// Fortran functions
extern "C" {
void GlobalMesh_create_mesh(int, int, int, const double[], const int[],
                            const int[], const int[], const int[], const int[],
                            const int[], const int[], const int[]);
void GlobalMesh_allocate_petrophysics();
void GlobalMesh_set_all_rocktypes();
void build_well_connectivity();
void GlobalMesh_mesh_bounding_box();
void GlobalMesh_compute_all_connectivies();
void GlobalMesh_set_frac();
void GlobalMesh_node_of_frac();
void GlobalMesh_frac_by_node();
void DefWell_make_compute_well_index();
void set_Peaceman_WI_threshold(double);
void unset_Peaceman_WI_threshold();
void GlobalMesh_allocate_rocktype();
void GlobalMesh_count_dirichlet_nodes();
}

namespace py = pybind11;

// FIXME: Retrieve id and coordinate types from meshtools (avoid double and int
// here)
void create_mesh(py::array_t<double, py::array::c_style> vertices,
                 py::array_t<int, py::array::c_style> cells_nodes_pointers,
                 py::array_t<int, py::array::c_style> cells_nodes_values,
                 py::array_t<int, py::array::c_style> cells_faces_pointers,
                 py::array_t<int, py::array::c_style> cells_faces_values,
                 py::array_t<int, py::array::c_style> faces_nodes_pointers,
                 py::array_t<int, py::array::c_style> faces_nodes_values) {
   std::vector<int> cellids;
   const std::size_t nb_cells = cells_nodes_pointers.shape(0) - 1;
   assert(nb_cells > 0);
   for (std::size_t i = 0; i != nb_cells; ++i) {
      cellids.emplace_back(0);
   }
   std::vector<int> faceids;
   const std::size_t nb_faces = faces_nodes_pointers.shape(0) - 1;
   assert(nb_faces > 0);
   for (std::size_t i = 0; i != nb_faces; ++i) {
      faceids.emplace_back(0);
   }
   GlobalMesh_create_mesh(
       vertices.shape(0), nb_cells, nb_faces, vertices.data(0, 0),
       cells_faces_pointers.data(0), cells_faces_values.data(0),
       cells_nodes_pointers.data(0), cells_nodes_values.data(0),
       faces_nodes_pointers.data(0), faces_nodes_values.data(0), cellids.data(),
       faceids.data());
}

// void create_mesh(
//	py::array_t<double, py::array::c_style> vertices,
//	py::tuple cells_nodes,
//	py::tuple cells_faces,
//	py::tuple faces_nodes
//)
//{
//	auto cells_nodes_pointers = cells_nodes[0];
//	auto cells_nodes_values = cells_nodes[1];
//	auto cells_faces_pointers = cells_faces[0];
//	auto cells_faces_values = cells_faces[1];
//	auto faces_nodes_pointers = faces_nodes[0];
//	auto faces_nodes_values = faces_nodes[1];
//	create_mesh(vertices,
//		cells_nodes_pointers, cells_nodes_values,
//		cells_faces_pointers, cells_faces_values,
//		faces_nodes_pointers, faces_nodes_values);
//}

void add_GlobalMesh_wrappers(py::module& module) {
   // module.def("create_mesh", [](
   //	py::array_t<double, py::array::c_style> vertices,
   //	py::tuple connectivity
   //	) {
   //	create_mesh(vertices, connectivity[0], connectivity[1],
   // connectivity[2]);
   //});
   module.def("create_mesh", &create_mesh);

   // This is only transitory
   // module.def("global_mesh_make_post_read", &GlobalMesh_make_post_read,
   // "Compute all well indices.");
   // module.def("global_mesh_make_post_read_fracture_and_dirBC",
   // &GlobalMesh_make_post_read_fracture_and_dirBC, "Set fractures and boudary
   // conditions.");
   module.def("global_mesh_allocate_petrophysics",
              &GlobalMesh_allocate_petrophysics,
              "Allocate porosity and permeability.");
   module.def("global_mesh_set_all_rocktypes", &GlobalMesh_set_all_rocktypes,
              "Set node rocktypes taking into account neighboring "
              "permeabilities/conductivities.");
   module.def("build_well_connectivity", &build_well_connectivity,
              "Compute well connectivity.");
   module.def("global_mesh_mesh_bounding_box", &GlobalMesh_mesh_bounding_box);
   module.def("global_mesh_compute_all_connectivies",
              &GlobalMesh_compute_all_connectivies);
   module.def("global_mesh_set_frac", &GlobalMesh_set_frac);
   module.def("global_mesh_node_of_frac", &GlobalMesh_node_of_frac);
   // module.def("global_mesh_set_dir_BC", &GlobalMesh_set_dir_BC);
   module.def("global_mesh_frac_by_node", &GlobalMesh_frac_by_node);
   module.def("global_mesh_allocate_rocktype", &GlobalMesh_allocate_rocktype,
              "Allocate NodeRocktype, FracRocktype, CellRocktype");
   module.def("global_mesh_count_dirichlet_nodes",
              &GlobalMesh_count_dirichlet_nodes);

   module.def("compute_well_indices", &DefWell_make_compute_well_index,
              "Compute all well indices.");
   module.def("set_peaceman_index_threshold", &set_Peaceman_WI_threshold);
   module.def("unset_peaceman_index_threshold", &unset_Peaceman_WI_threshold);
}
