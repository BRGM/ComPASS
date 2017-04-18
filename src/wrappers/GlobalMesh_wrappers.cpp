// Fortran functions
extern "C"
{
	void GlobalMesh_build_cartesian_grid(double, double, double, double, double, double, int, int, int);
	void GlobalMesh_make_post_read();
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

	module.def("compute_well_indices", &DefWell_make_compute_well_index, "Compute all well indices.");

}
