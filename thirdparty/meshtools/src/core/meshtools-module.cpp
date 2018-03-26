#include "meshtools-wrapper.h"

PYBIND11_MODULE(_MeshTools, module)
{

	module.doc() = "pybind11 generic mesh tools (quick and dirty!!!)";
	add_mesh_tools(module);

}
