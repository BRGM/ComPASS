// Fortran functions
extern "C"
{
	void GlobalMesh_build_cartesian_grid(double, double, double, double, double, double, int, int, int);
	void GlobalMesh_create_mesh(int, int, int, const double[], const int[], const int[], const int[], const int[], const int[], const int[], const int[], const int[]);
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
	void GlobalMesh_allocate_id_nodes();
	void GlobalMesh_count_dirichlet_nodes();
	void GlobalMesh_set_cartesian_mesh();
	void GlobalMesh_set_hexahedron_mesh();
	void GlobalMesh_set_tetrahedron_mesh();
	void GlobalMesh_set_wedge_mesh();
}

#include "GlobalMesh_wrappers.h"
#include <pybind11/numpy.h>

// FIXME: Retrive id and coordinate types from meshtools (avoid double and int here)
void create_mesh(
	py::array_t<double, py::array::c_style> vertices,
	py::array_t<int, py::array::c_style> cells_nodes_pointers,
	py::array_t<int, py::array::c_style> cells_nodes_values,
	py::array_t<int, py::array::c_style> cells_faces_pointers,
	py::array_t<int, py::array::c_style> cells_faces_values,
	py::array_t<int, py::array::c_style> faces_nodes_pointers,
	py::array_t<int, py::array::c_style> faces_nodes_values
)
{
	std::vector<int> cellids;
	const std::size_t nb_cells = cells_nodes_pointers.shape(0) - 1;
	assert(nb_cells > 1);
	for (std::size_t i=0; i != nb_cells; ++i) {
		cellids.emplace_back(0);
	}
	std::vector<int> faceids;
	const std::size_t nb_faces = faces_nodes_pointers.shape(0) - 1;
	assert(nb_faces > 0);
	for (std::size_t i = 0; i != nb_faces; ++i) {
		faceids.emplace_back(0);
	}
	GlobalMesh_create_mesh(
		vertices.shape(0), nb_cells, nb_faces,
		vertices.data(0, 0),
		cells_faces_pointers.data(0), cells_faces_values.data(0),
		cells_nodes_pointers.data(0), cells_nodes_values.data(0),
		faces_nodes_pointers.data(0), faces_nodes_values.data(0),
		cellids.data(), faceids.data()
	);
}

//void create_mesh(
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
		"Build a cartesian grid. This routine must be called by the master process.");
	
	//module.def("create_mesh", [](
	//	py::array_t<double, py::array::c_style> vertices,
	//	py::tuple connectivity
	//	) {
	//	create_mesh(vertices, connectivity[0], connectivity[1], connectivity[2]);
	//});
	module.def("create_mesh", &create_mesh);

	// This is only transitory
	module.def("global_mesh_make_post_read", &GlobalMesh_make_post_read, "Compute all well indices.");
	module.def("global_mesh_make_post_read_fracture_and_dirBC", &GlobalMesh_make_post_read_fracture_and_dirBC, "Set fractures and boudary conditions.");
	module.def("global_mesh_make_post_read_set_poroperm", &GlobalMesh_make_post_read_set_poroperm, "Set porosity and permeability.");
	module.def("global_mesh_make_post_read_well_connectivity_and_ip", &GlobalMesh_make_post_read_well_connectivity_and_ip, "Compute well connectivity and PI.");
	module.def("global_mesh_mesh_bounding_box", &GlobalMesh_mesh_bounding_box);
	module.def("global_mesh_compute_all_connectivies", &GlobalMesh_compute_all_connectivies);
	module.def("global_mesh_set_frac", &GlobalMesh_set_frac);
	module.def("global_mesh_node_of_frac", &GlobalMesh_node_of_frac);
	//module.def("global_mesh_set_dir_BC", &GlobalMesh_set_dir_BC);
	module.def("global_mesh_frac_by_node", &GlobalMesh_frac_by_node);
	module.def("global_mesh_allocate_id_nodes", &GlobalMesh_allocate_id_nodes);
	module.def("global_mesh_count_dirichlet_nodes", &GlobalMesh_count_dirichlet_nodes);
	module.def("global_mesh_set_cartesian_mesh", &GlobalMesh_set_cartesian_mesh);
	module.def("global_mesh_set_hexahedron_mesh", &GlobalMesh_set_hexahedron_mesh);
	module.def("global_mesh_set_tetrahedron_mesh", &GlobalMesh_set_tetrahedron_mesh);
	module.def("global_mesh_set_wedge_mesh", &GlobalMesh_set_wedge_mesh);

	module.def("compute_well_indices", &DefWell_make_compute_well_index, "Compute all well indices.");

}