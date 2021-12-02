//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>

// FIXME: should be defined with common types/defines
//        in ComPASS namespace
typedef int ComPASS_Fortran_size_type;

// Fortran functions
extern "C" {
double get_gravity();
void set_gravity(double);
#ifdef _WITH_FREEFLOW_STRUCTURES_
double get_atm_pressure();
void set_atm_pressure(double);
double get_atm_temperature();
void set_atm_temperature(double);
double get_atm_flux_radiation();
void set_atm_flux_radiation(double);
double get_soil_emissivity();
void set_soil_emissivity(double);
double get_atm_rain_flux();
void set_atm_rain_flux(double);
#endif
double get_rock_volumetric_heat_capacity();
void set_rock_volumetric_heat_capacity(double);
double get_fracture_thickness();
void set_fracture_thickness(double);
}

namespace py = pybind11;

void add_global_variables_wrappers(py::module& module) {
   module.def("get_gravity", &get_gravity);
   module.def("set_gravity", &set_gravity);
#ifdef _WITH_FREEFLOW_STRUCTURES_
   module.def("get_atm_pressure", &get_atm_pressure);
   module.def("set_atm_pressure", &set_atm_pressure);
   module.def("get_atm_temperature", &get_atm_temperature);
   module.def("set_atm_temperature", &set_atm_temperature);
   module.def("get_atm_flux_radiation", &get_atm_flux_radiation);
   module.def("set_atm_flux_radiation", &set_atm_flux_radiation);
   module.def("get_soil_emissivity", &get_soil_emissivity);
   module.def("set_soil_emissivity", &set_soil_emissivity);
   module.def("get_atm_rain_flux", &get_atm_rain_flux);
   module.def("set_atm_rain_flux", &set_atm_rain_flux);
#endif
   module.def("get_rock_volumetric_heat_capacity",
              &get_rock_volumetric_heat_capacity);
   module.def("set_rock_volumetric_heat_capacity",
              &set_rock_volumetric_heat_capacity);
   module.def("get_fracture_thickness", &get_fracture_thickness);
   module.def("set_fracture_thickness", &set_fracture_thickness);
   module.def("has_energy_transfer_enabled", []() {
#ifdef _THERMIQUE_
      return true;
#else
        return false;
#endif
   });
}
