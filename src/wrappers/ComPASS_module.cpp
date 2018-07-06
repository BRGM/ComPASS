//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

#include "preprocessor_utils.h"
#include "PyBuffer_wrappers.h"
#include "NN_wrappers.h"
#include "COC_wrappers.h"
#include "GlobalMesh_wrappers.h"
#include "GlobalVariables_wrappers.h"
#include "IncCV_wrappers.h"
#include "MeshUtilities_wrappers.h"
#include "TimeLoop_wrappers.h"
#include "DebugUtils_wrappers.h"
#include "Model_wrappers.h"
#include "Flux_wrappers.h"
#include "Well_wrappers.h"

PYBIND11_MODULE(ComPASS_CONFIGURATION_NAME, module)
{

	module.doc() = "pybind11 ComPASS library interface for " TOSTRING(ComPASS_CONFIGURATION_NAME) " configuration";

	add_pybuffer_wrappers(module);
	add_NN_wrappers(module);
	add_coc_wrappers(module);
	add_GlobalMesh_wrappers(module);
	add_mesh_utilities_wrappers(module);
	add_well_wrappers(module);
	add_global_variables_wrappers(module);
    add_IncCV_wrappers(module);
    add_Model_wrappers(module);
    add_flux_wrappers(module);
    add_time_loop_wrappers(module);
    add_debug_utils_wrappers(module);

}
