//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

namespace py = pybind11;

void add_pybuffer_wrappers(py::module &);
void add_NN_wrappers(py::module &);
void add_coc_wrappers(py::module &);
void add_GlobalMesh_wrappers(py::module &);
void add_mesh_utilities_wrappers(py::module &);
void add_mesh_schema_wrappers(py::module &);
void add_well_wrappers(py::module &);
void add_global_variables_wrappers(py::module &);
void add_IncCV_wrappers(py::module &);
void add_flux_wrappers(py::module &);
void add_time_loop_wrappers(py::module &);
void add_debug_utils_wrappers(py::module &);
void add_model_wrappers(py::module &);
void add_Residu_wrappers(py::module &);
void add_Metis_wrapper(py::module &);
void add_SolvePetsc_wrappers(py::module &);
void add_LinearSystem_wrapper(py::module &);
void add_Jacobian_wrappers(py::module &);
void add_SyncPetsc_wrappers(py::module &);
void add_VAGFrac_wrappers(py::module &);
void add_freeflow_wrappers(py::module &);
void add_mswell_wrappers(py::module &module);
void add_ResiduMSWells_wrappers(py::module &module);
void add_LinearSystemMSWells_wrapper(py::module &module);
void add_physical_properties_wrappers(py::module &module);

// FIXME: this should be elsewhere as it internals and not wrappers
void add_simulation_pyinternals(py::module &);
void add_petrophysics_pyinternals(py::module &);

#include "preprocessor_utils.h"

PYBIND11_MODULE(ComPASS_CONFIGURATION_NAME, module) {
   module.doc() = "pybind11 ComPASS library interface for " TOSTRING(
       ComPASS_CONFIGURATION_NAME) " configuration";

   add_pybuffer_wrappers(module);
   add_NN_wrappers(module);
   add_coc_wrappers(module);
   add_GlobalMesh_wrappers(module);
   add_mesh_utilities_wrappers(module);
   add_mesh_schema_wrappers(module);
   add_well_wrappers(module);
   add_global_variables_wrappers(module);
   add_IncCV_wrappers(module);
   add_flux_wrappers(module);
   add_time_loop_wrappers(module);
   add_debug_utils_wrappers(module);
   add_model_wrappers(module);
   add_Residu_wrappers(module);
   add_Metis_wrapper(module);
   add_SolvePetsc_wrappers(module);
   add_LinearSystem_wrapper(module);
   add_Jacobian_wrappers(module);
   add_SyncPetsc_wrappers(module);
   add_VAGFrac_wrappers(module);
   add_freeflow_wrappers(module);
   add_mswell_wrappers(module);
   add_ResiduMSWells_wrappers(module);
   add_LinearSystemMSWells_wrapper(module);
   add_physical_properties_wrappers(module);

   // FIXME: this should be elsewhere as it is more an internal than a wrapper
   add_simulation_pyinternals(module);
   add_petrophysics_pyinternals(module);
}
