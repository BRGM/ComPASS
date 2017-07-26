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

#include "meshtools.h"
#include "meshtools-wrapper.h"

namespace MT = MeshTools;

template <typename Mesh>
auto create_mesh(const Mesh& mesh)
{
	const auto& vertices = mesh.vertices;
	const auto& cells = mesh.connectivity.cells;
	const auto& faces = mesh.connectivity.faces;
	const auto cellnodes = MT::FSCoC_as_COC(cells.nodes);
	auto cellnodes_pointers = std::get<0>(cellnodes);
	auto cellnodes_values = std::get<1>(cellnodes);
	const auto cellfaces = MT::FSCoC_as_COC(cells.faces);
	auto cellfaces_pointers = std::get<0>(cellfaces);
	auto cellfaces_values = std::get<1>(cellfaces);
	const auto facenodes = MT::FSCoC_as_COC(faces.nodes);
	auto facenodes_pointers = std::get<0>(facenodes);
	auto facenodes_values = std::get<1>(facenodes);
	std::vector<int> cellids;
	std::size_t n = cells.nb();
	for (; n != 0; --n) {
		cellids.emplace_back(0);
	}
	std::vector<int> faceids;
	n = faces.nb();
	for (; n != 0; --n) {
		faceids.emplace_back(0);
	}
	GlobalMesh_create_mesh(
		vertices.size(), cells.nb(), faces.nb(),
		vertices.data()->data(),
		cellfaces_pointers.data(), cellfaces_values,
		cellnodes_pointers.data(), cellnodes_values,
		facenodes_pointers.data(), facenodes_values,
		cellids.data(), faceids.data()
	);
}

void add_GlobalMesh_wrappers(py::module& module)
{

	// add meshtools submodule
	auto mesh_tools_module = module.def_submodule("MeshTools", "MeshTools submodules.");
	add_mesh_tools(mesh_tools_module);

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

	// CHECKME: The weird conversion is due to a gcc bug
	// cf. https://stackoverflow.com/questions/45077622/using-a-template-function-pointer-inside-another-template-function?noredirect=1#comment77158283_45077622
	module.def("create_mesh", (decltype(&create_mesh<MT::TetMesh>))(&create_mesh<MT::TetMesh>), "Creates a tet mesh.");
	module.def("create_mesh", (decltype(&create_mesh<MT::HexMesh>))(&create_mesh<MT::HexMesh>), "Creates a hex mesh.");
	
	//void GlobalMesh_create_mesh(int, int, int, double[], int[], int[], int[], int[], int[], int[], int[], int[]);
	//subroutine GlobalMesh_create_mesh_from_C(nbnodes, nbcells, nbfaces, &
	//	nodes, &
	//	cell_faces_ptr, cell_faces_val, &
	//	cell_nodes_ptr, cell_nodes_val, &
	//	face_nodes_ptr, face_nodes_val, &
	//	cell_id, face_id) &
	//	bind(C, name = "GlobalMesh_create_mesh")

	//	integer(c_int), value, intent(in) ::nbnodes
	//	integer(c_int), value, intent(in) ::nbcells
	//	integer(c_int), value, intent(in) ::nbfaces
	//	real(c_double), dimension(3, nbnodes), intent(in) ::nodes
	//	integer(c_int), dimension(nbcells + 1), intent(in) ::cell_faces_ptr
	//	integer(c_int), dimension(cell_faces_ptr(nbcells + 1)), intent(in) ::cell_faces_val
	//	integer(c_int), dimension(nbcells + 1), intent(in) ::cell_nodes_ptr
	//	integer(c_int), dimension(cell_nodes_ptr(nbcells + 1)), intent(in) ::cell_nodes_val
	//	integer(c_int), dimension(nbfaces + 1), intent(in) ::face_nodes_ptr
	//	integer(c_int), dimension(face_nodes_ptr(nbfaces + 1)), intent(in) ::face_nodes_val
	//	integer(c_int), dimension(nbcells), intent(in) ::cell_id
	//	integer(c_int), dimension(nbfaces), intent(in) ::face_id

	//	call GlobalMesh_create_mesh(nodes, &
	//		cell_faces_ptr, cell_faces_val, &
	//		cell_nodes_ptr, cell_nodes_val, &
	//		face_nodes_ptr, face_nodes_val, &
	//		cell_id, face_id, &
	//		.true.)


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
