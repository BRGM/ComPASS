#include "ArrayWrapper.h"
#include "PyBuffer_wrappers.h"
#include "MeshUtilities.h"

// Fortran functions
extern "C"
{
	void retrieve_vertices(ArrayWrapper&);
	void retrieve_mesh_connectivity(MeshConnectivity&);
	void retrieve_id_faces(ArrayWrapper&);
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

}