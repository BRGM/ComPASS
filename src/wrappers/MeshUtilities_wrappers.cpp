#include "ArrayWrapper.h"
#include "XArrayWrapper.h"
#include "PyXArrayWrapper.h"
#include "PyArrayWrapper.h"
#include "PyBuffer_wrappers.h"
#include "MeshUtilities.h"


// FUTURE: Code information bitwise?
struct NodeInfo
{
	char proc; // 'o' / 'g': own / ghost
	char frac; // 'y' / 'n' : node in fracture / not in fracture
			   // FIXME: Name is to be changed
	char pressure; // 'd' / 'n' / 'i' : dirichlet / neumann / interior for the pressure
				   // FIXME: preprocessor directives to be removed!
#ifdef _THERMIQUE_
				   // FIXME: Name is to be changed
	char temperature; // 'd' / 'n' / 'i' : dirichlet / neumann / interior for the temperature
#endif
};

struct Point
{
	double x, y, z;
};

// Fortran functions
extern "C"
{
	void retrieve_vertices(XArrayWrapper<Point>&);
	void retrieve_nodeflags(XArrayWrapper<int>&);
	void retrieve_global_nodeflags(XArrayWrapper<int>&);
	void retrieve_cellflags(XArrayWrapper<int>&);
	void retrieve_global_cellflags(XArrayWrapper<int>&);
	void retrieve_celltypes(XArrayWrapper<int8_t>&);
	void retrieve_global_celltypes(XArrayWrapper<int8_t>&);
	void retrieve_facetypes(XArrayWrapper<int8_t>&);
	void retrieve_global_facetypes(XArrayWrapper<int8_t>&);
	void retrieve_faceflags(XArrayWrapper<int>&);
	void retrieve_global_faceflags(XArrayWrapper<int>&);
	void retrieve_global_vertices(XArrayWrapper<Point>&);
	void retrieve_global_mesh_connectivity(MeshConnectivity&);
	void retrieve_mesh_connectivity(MeshConnectivity&);
	void retrieve_id_faces(ArrayWrapper&);
	void retrieve_cell_porosity(ArrayWrapper&);
	void retrieve_face_porosity(ArrayWrapper&);
	void retrieve_cell_permeability(ArrayWrapper&);
	void retrieve_face_permeability(ArrayWrapper&);
	void retrieve_global_id_node(XArrayWrapper<NodeInfo>&);
	void retrieve_id_node(XArrayWrapper<NodeInfo>&);
	void retrieve_frac_face_id(XArrayWrapper<int>&);
	void retrieve_face_frac_id(XArrayWrapper<int>&);
	void retrieve_nb_cells_own(XArrayWrapper<int>&);
	void retrieve_nb_faces_own(XArrayWrapper<int>&);
	void retrieve_nb_nodes_own(XArrayWrapper<int>&);
	void retrieve_nb_fractures_own(XArrayWrapper<int>&);
}

#include "MeshUtilities_wrappers.h"

void add_mesh_utilities_wrappers(py::module& module)
{

	PYBIND11_NUMPY_DTYPE(Point, x, y, z);
	add_array_wrapper(module, "global_vertices", retrieve_global_vertices);
	add_array_wrapper(module, "global_nodeflags", retrieve_global_nodeflags);
	add_array_wrapper(module, "global_cellflags", retrieve_global_cellflags);
	add_array_wrapper(module, "global_faceflags", retrieve_global_faceflags);
	add_array_wrapper(module, "global_celltypes", retrieve_global_celltypes);
	add_array_wrapper(module, "global_facetypes", retrieve_global_facetypes);
	add_array_wrapper(module, "vertices", retrieve_vertices);
	add_array_wrapper(module, "nodeflags", retrieve_nodeflags);
	add_array_wrapper(module, "cellflags", retrieve_cellflags);
	add_array_wrapper(module, "faceflags", retrieve_faceflags);
	add_array_wrapper(module, "celltypes", retrieve_celltypes);
	add_array_wrapper(module, "facetypes", retrieve_facetypes);
	add_array_wrapper(module, "face_frac_id", retrieve_face_frac_id);
	add_array_wrapper(module, "frac_face_id", retrieve_frac_face_id);
	add_array_wrapper(module, "nb_cells_own", retrieve_nb_cells_own);
	add_array_wrapper(module, "nb_faces_own", retrieve_nb_faces_own);
	add_array_wrapper(module, "nb_nodes_own", retrieve_nb_nodes_own);
	add_array_wrapper(module, "nb_fractures_own", retrieve_nb_fractures_own);

	module.def("get_id_faces_buffer",
		[]() { return retrieve_buffer<IntBuffer>(retrieve_id_faces); },
		"Get faces integer flag. Can be used to specify fracture faces setting the flag to -2.");

	module.def("get_cell_porosity_buffer",
		[]() { return retrieve_buffer<DoubleBuffer>(retrieve_cell_porosity); }
	);

	module.def("get_face_porosity_buffer",
		[]() { return retrieve_buffer<DoubleBuffer>(retrieve_face_porosity); }
	);

	module.def("get_cell_permeability_buffer",
		[]() { return retrieve_buffer<TensorBuffer>(retrieve_cell_permeability); }
	);

	module.def("get_face_permeability_buffer",
		[]() { return retrieve_buffer<DoubleBuffer>(retrieve_face_permeability); }
	);

	py::class_<MeshConnectivity>(module, "MeshConnectivity")
		.def_readwrite("NodebyCell", &MeshConnectivity::NodebyCell)
		.def_readwrite("NodebyFace", &MeshConnectivity::NodebyFace)
		.def_readwrite("FacebyCell", &MeshConnectivity::FacebyCell)
		.def_readwrite("CellbyNode", &MeshConnectivity::CellbyNode)
		.def_readwrite("CellbyFace", &MeshConnectivity::CellbyFace)
		.def_readwrite("CellbyCell", &MeshConnectivity::CellbyCell);

	module.def("get_global_connectivity", []() {
		MeshConnectivity connectivity;
		retrieve_global_mesh_connectivity(connectivity);
		return connectivity;
	},
		"Get global mesh connectivity."
		);

	module.def("get_connectivity", []() {
		MeshConnectivity connectivity;
		retrieve_mesh_connectivity(connectivity);
		return connectivity;
	},
		"Get local mesh connectivity."
		);

	// FIXME: preprocessor directives to be removed!
#ifdef _THERMIQUE_
	PYBIND11_NUMPY_DTYPE(NodeInfo, proc, frac, pressure, temperature);
	#else
	PYBIND11_NUMPY_DTYPE(NodeInfo, proc, frac, pressure);
#endif

	add_array_wrapper(module, "global_node_info", retrieve_global_id_node);
	add_array_wrapper(module, "node_info", retrieve_id_node);

}