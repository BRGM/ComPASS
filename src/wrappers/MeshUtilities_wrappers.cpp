#include "MeshUtilities.h"

// Fortran functions
extern "C"
{
	void retrieve_vertices(Vertices&);
	void retrieve_mesh_connectivity(MeshConnectivity&);
}

#include "MeshUtilities_wrappers.h"

void add_mesh_utilities_wrappers(py::module& module)
{

	py::class_<Vertices>(module, "VerticesHandle", py::buffer_protocol())
		.def_buffer([](Vertices &V) -> py::buffer_info {
		return py::buffer_info(
			V.p,                               /* Pointer to buffer */
			sizeof(double),                          /* Size of one scalar */
			py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
			2,                                      /* Number of dimensions */
			{ V.nb_points, 3 },                 /* Buffer dimensions */
			{ sizeof(double) * 3,             /* Strides (in bytes) for each index */
			sizeof(double) }
		);
	});

	module.def("get_vertices", []() {
		Vertices V;
		retrieve_vertices(V);
		return V;
	},
		"Get node coordinates."
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

}