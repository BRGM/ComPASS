//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <cstddef>
#include <limits>

#include "CTVector.h"
#include "NewtonIncrements.h"
#include "StringWrapper.h"

// Fortran functions
extern "C" {
void IncCV_SaveIncPreviousTimeStep();
void IncCV_LoadIncPreviousTimeStep();
void IncCVWells_estimate_producers_density(const bool);
void IncCVWells_UpdatePressureDrop();
void IncCVWells_UpdateWellPressures();
double Newton_compute_relaxation(
    const NewtonIncrements::Pointers<const double>);
void IncCV_NewtonIncrement(const NewtonIncrements::Pointers<const double>,
                           const double);
void DirichletContribution_update();
void IncPrimSecd_compute();
void IncPrimSecdFreeFlow_compute();
void LoisThermoHydro_compute_phase_pressures();
void LoisThermoHydro_compute_phase_pressures_derivatives();
void LoisThermoHydro_compute();
void Flux_DarcyFlux_Cell();
void Flux_DarcyFlux_Frac();
void Flux_FourierFlux_Cell();
void Flux_FourierFlux_Frac();
void Residu_reset_history();
void Residu_compute(double);
double Residu_RelativeNorm_local_closure();
void IncPrimSecd_PrimToSecd(NewtonIncrements::Pointers<double>);
void NN_flash_all_control_volumes();
void DefFlashWells_TimeFlash_producers_single_phase();
void DefFlashWells_TimeFlash_producers_two_phases();
void DefFlashWells_TimeFlash_injectors();
void pass_and_dump_array(double*, std::size_t*);
}

namespace py = pybind11;

void add_time_loop_wrappers(py::module& module) {
   module.def("IncCV_SaveIncPreviousTimeStep", &IncCV_SaveIncPreviousTimeStep);
   module.def("IncCV_LoadIncPreviousTimeStep", &IncCV_LoadIncPreviousTimeStep);
   module.def("IncCVWells_estimate_producers_density",
              &IncCVWells_estimate_producers_density);
   module.def("IncCVWells_UpdatePressureDrop", &IncCVWells_UpdatePressureDrop);
   module.def("IncCVWells_UpdateWellPressures",
              &IncCVWells_UpdateWellPressures);
   module.def("DirichletContribution_update", &DirichletContribution_update);
   module.def("IncPrimSecd_update_secondary_dependencies", []() {
      IncPrimSecd_compute();
#ifdef _WITH_FREEFLOW_STRUCTURES_
      IncPrimSecdFreeFlow_compute();
#endif  // _WITH_FREEFLOW_STRUCTURES_
   });
   module.def("LoisThermoHydro_compute_phase_pressures",
              &LoisThermoHydro_compute_phase_pressures);
   module.def("LoisThermoHydro_compute_phase_pressures_derivatives",
              &LoisThermoHydro_compute_phase_pressures_derivatives);
   module.def("LoisThermoHydro_compute", &LoisThermoHydro_compute);
   module.def("Flux_DarcyFlux_Cell", &Flux_DarcyFlux_Cell);
   module.def("Flux_DarcyFlux_Frac", &Flux_DarcyFlux_Frac);
   module.def("Flux_FourierFlux_Cell", &Flux_FourierFlux_Cell);
   module.def("Flux_FourierFlux_Frac", &Flux_FourierFlux_Frac);
   module.def("Residu_compute", &Residu_compute);
   module.def("Residu_reset_history", &Residu_reset_history);
   module.def("Residu_RelativeNorm_local_closure",
              &Residu_RelativeNorm_local_closure);
   module.def("NN_flash_all_control_volumes", &NN_flash_all_control_volumes);
   module.def("DefFlashWells_TimeFlash_producers_single_phase",
              &DefFlashWells_TimeFlash_producers_single_phase);
   module.def("DefFlashWells_TimeFlash_producers_two_phases",
              &DefFlashWells_TimeFlash_producers_two_phases);
   module.def("DefFlashWells_TimeFlash_injectors",
              &DefFlashWells_TimeFlash_injectors);

   auto as_primary_array = [](std::vector<double>& v) {
      assert(v.size() < std::numeric_limits<py::ssize_t>::max());
      auto res = py::array_t<double, py::array::c_style>{
          static_cast<py::ssize_t>(v.size()), v.data()};
      res.attr("shape") = py::make_tuple(-1, NewtonIncrements::npv());
      return res;
   };

   py::class_<NewtonIncrements>(module, "NewtonIncrements")
       .def(py::init())
       .def("init", &NewtonIncrements::init)
       .def(
           "nodes",
           [as_primary_array](NewtonIncrements& self) {
              return as_primary_array(self.nodes);
           },
           py::keep_alive<0, 1>())
       .def(
           "fractures",
           [as_primary_array](NewtonIncrements& self) {
              return as_primary_array(self.fractures);
           },
           py::keep_alive<0, 1>())
       .def(
           "cells",
           [as_primary_array](NewtonIncrements& self) {
              return as_primary_array(self.cells);
           },
           py::keep_alive<0, 1>());

   py::class_<CTVector>(module, "CTVector")
       .def(py::init())
       .def(
           "as_array",
           [](CTVector& self) {
              return py::array_t<double, py::array::c_style>{
                  static_cast<py::ssize_t>(self.npv), self.values};
           },
           py::keep_alive<0, 1>());

   module.def("Newton_compute_relaxation",
              [](const NewtonIncrements& increments) {
                 return Newton_compute_relaxation(increments.pointers());
              });

   module.def("IncCV_NewtonIncrement",
              [](const NewtonIncrements& increments, const double relaxation) {
                 IncCV_NewtonIncrement(increments.pointers(), relaxation);
              });

   module.def("IncPrimSecd_PrimToSecd", [](NewtonIncrements& increments) {
      //     py::print("LT P2S sizes:", nb_primary_variables(), nb_nodes(),
      //			    nb_fractures(), nb_cells(), nb_injectors(),
      // nb_producers());
      IncPrimSecd_PrimToSecd(increments.pointers());
   });
}
