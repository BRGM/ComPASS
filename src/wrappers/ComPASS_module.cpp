#include <pybind11/pybind11.h>

#include "NN_wrappers.h"
#include "COC_wrappers.h"
#include "GlobalMesh_wrappers.h"
#include "MeshUtilities_wrappers.h"
#include "Well_wrappers.h"

PYBIND11_PLUGIN(ComPASS)
{

	py::module module("ComPASS", "pybind11 ComPASS library interface");
	
	add_NN_wrappers(module);
	add_coc_wrappers(module);
	add_GlobalMesh_wrappers(module);
	add_mesh_utilities_wrappers(module);
	add_well_wrappers(module);

	return module.ptr();

}
