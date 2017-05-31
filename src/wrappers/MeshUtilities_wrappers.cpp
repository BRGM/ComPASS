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

// Fortran functions
extern "C"
{
	void retrieve_vertices(ArrayWrapper&);
	void retrieve_mesh_connectivity(MeshConnectivity&);
	void retrieve_id_faces(ArrayWrapper&);
	void retrieve_cell_porosity(ArrayWrapper&);
	void retrieve_face_porosity(ArrayWrapper&);
	void retrieve_cell_permeability(ArrayWrapper&);
	void retrieve_face_permeability(ArrayWrapper&);
	void retrieve_global_id_node(XArrayWrapper<NodeInfo>&);
	void retrieve_id_node(XArrayWrapper<NodeInfo>&);
}

#include "MeshUtilities_wrappers.h"

void add_mesh_utilities_wrappers(py::module& module)
{

	module.def("get_vertices_buffer",
		[]() { return retrieve_buffer<CoordinatesBuffer>(retrieve_vertices); },
		"Get node coordinates.");

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

	module.def("get_connectivity", []() {
		MeshConnectivity connectivity;
		retrieve_mesh_connectivity(connectivity);
		//std::cout << connectivity << std::endl;
		return connectivity;
	},
		"Get mesh connectivity."
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