// Fortran functions
extern "C"
{
  void GlobalMesh_build_cartesian_grid(double, double, double, double, double, double, int, int, int);
  void GlobalMesh_make_post_read();
  void GlobalMesh_make_post_read_fracture_and_dirBC();
  void GlobalMesh_make_post_read_set_poroperm();
  void GlobalMesh_make_post_read_well_connectivity_and_ip();
  void GlobalMesh_mesh_bounding_box();
  void GlobalMesh_compute_all_connectivies();
  void GlobalMesh_set_frac();
  void GlobalMesh_node_of_frac();
  void GlobalMesh_set_dir_BC();
  void GlobalMesh_frac_by_node();
  void DefWell_make_compute_well_index();
}

#include "GlobalMesh_wrappers.h"

void add_GlobalMesh_wrappers(py::module& module)
{

	module.def("build_grid",
		[](py::object shape, py::object extent, py::object origin) {
		if (origin.is_none()) origin = py::make_tuple(0, 0, 0);
		if (extent.is_none()) extent = py::make_tuple(1, 1, 1);
		auto shape_tuple = shape.cast<py::tuple>();
		auto extent_tuple = extent.cast<py::tuple>();
		auto origin_tuple = origin.cast<py::tuple>();
		GlobalMesh_build_cartesian_grid(
			origin_tuple[0].cast<double>(), origin_tuple[1].cast<double>(), origin_tuple[2].cast<double>(),
			extent_tuple[0].cast<double>(), extent_tuple[1].cast<double>(), extent_tuple[2].cast<double>(),
			shape_tuple[0].cast<int>(), shape_tuple[1].cast<int>(), shape_tuple[2].cast<int>()
			);
	},
		py::arg("shape"), py::arg("extent") = py::none{}, py::arg("origin") = py::none{},
		"Build a cartesian grid. This routine must be called by the master process." );

	// This is only transitory
	module.def("global_mesh_make_post_read", &GlobalMesh_make_post_read, "Compute all well indices.");
	module.def("global_mesh_make_post_read_fracture_and_dirBC", &GlobalMesh_make_post_read_fracture_and_dirBC, "Set fractures and boudary conditions.");
	module.def("global_mesh_make_post_read_set_poroperm", &GlobalMesh_make_post_read_set_poroperm, "Set porosity and permeability.");
	module.def("global_mesh_make_post_read_well_connectivity_and_ip", &GlobalMesh_make_post_read_well_connectivity_and_ip, "Compute well connectivity and PI.");
	module.def("global_mesh_mesh_bounding_box", &GlobalMesh_mesh_bounding_box);
	module.def("global_mesh_compute_all_connectivies", &GlobalMesh_compute_all_connectivies);
	module.def("global_mesh_set_frac", &GlobalMesh_set_frac);
	module.def("global_mesh_node_of_frac", &GlobalMesh_node_of_frac);
	module.def("global_mesh_set_dir_BC", &GlobalMesh_set_dir_BC);
	module.def("global_mesh_frac_by_node", &GlobalMesh_frac_by_node);

	module.def("compute_well_indices", &DefWell_make_compute_well_index, "Compute all well indices.");

}
