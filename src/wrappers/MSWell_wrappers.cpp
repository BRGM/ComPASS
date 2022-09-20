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
#include <functional>
#include <iostream>
#include <set>
#include <sstream>

namespace py = pybind11;

// Fortran functions
#include "MeshSchema.fh"
#include "NewtonIncrements.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;

/** Common data structure shared by injectors and producers. */
struct Fortran_mswell_data {
   typedef std::array<double, NC> Component_vector;
   int id;
   double radius;
   double maximum_pressure;
   double minimum_pressure;
   double imposed_flowrate;
   Component_vector injection_composition;
   double injection_temperature;
   char operating_code;
   char well_type;
};

extern "C" {

void IncCVMSWells_copy_states_from_reservoir();
void IncPrimSecdMSWells_compute();
void MSWellsData_init();
void LoisThermoHydroMSWells_compute();
void VSHydroMSWells_init();
void VSHydroMSWells_compute();
void ResiduMSWells_compute(double);
void ResiduMSWells_reset_history();
void LeafMSWells_init_data();
void IncPrimSecdMSWells_PrimToSecd(NewtonIncrements::Pointers<double>);
double Newton_compute_relaxation_mswells(
    const NewtonIncrements::Pointers<const double>);
double Newton_get_last_max_inc_mswells();
// void IncCVMSWells_print_info_to_file();
void JacobianMSWells_ComputeJacSm(double, bool);
void JacobianMSWells_print_LA_info_to_file(double, int);
void JacobianMSWells_print_IP_info_to_file(double);
void IncCVMSWells_print_info_to_file();
void IncCVMSWells_set_compute_coupling(bool);

Fortran_mswell_data* get_mswells_data();
std::size_t nb_mswells();
std::size_t number_of_own_mswells();
}

struct MSWell {};

void add_mswell_wrappers(py::module& module) {
   py::class_<MSWell>(module, "MSWell").def(py::init<>());

   module.def("mswells_copy_states_from_reservoir",
              &IncCVMSWells_copy_states_from_reservoir,
              "For each vertex of all  multi-segmented producer wells, copy "
              "the  Coats variables from the reservoir");

   module.def("IncCVMSWells_set_compute_coupling",
              &IncCVMSWells_set_compute_coupling,
              "Compute coupling between resevoir and mswells");

   module.def("mswells_init_edge_data", &MSWellsData_init);
   module.def("IncPrimSecdMSWells_compute", &IncPrimSecdMSWells_compute);
   module.def("LoisThermoHydroMSWells_compute",
              &LoisThermoHydroMSWells_compute);

   module.def("mswells_init_leaf_data", &LeafMSWells_init_data);

   module.def("VSHydroMSWells_init", &VSHydroMSWells_init);
   module.def("VSHydroMSWells_compute", &VSHydroMSWells_compute);
   module.def("ResiduMSWells_compute", &ResiduMSWells_compute);
   module.def("ResiduMSWells_reset_history", &ResiduMSWells_reset_history);

   module.def("IncPrimSecdMSWells_PrimToSecd",
              [](NewtonIncrements& increments) {
                 IncPrimSecdMSWells_PrimToSecd(increments.pointers());
              });

   module.def(
       "Newton_compute_relaxation_mswells",
       [](const NewtonIncrements& increments) {
          return Newton_compute_relaxation_mswells(increments.pointers());
       });

   module.def("Newton_get_last_max_inc_mswells",
              []() { return Newton_get_last_max_inc_mswells(); });

   module.def("IncCVMSWells_print_info_to_file",
              &IncCVMSWells_print_info_to_file);

   module.def("JacobianMSWells_ComputeJacSm", &JacobianMSWells_ComputeJacSm);
   module.def("JacobianMSWells_print_LA_info_to_file",
              &JacobianMSWells_print_LA_info_to_file);

   module.def("JacobianMSWells_print_IP_info_to_file",
              &JacobianMSWells_print_IP_info_to_file);
   /////////////////////////////////////////////////////////////////////////
   py::class_<Fortran_mswell_data>(module, "MSWellData")
       .def_readonly("id", &Fortran_mswell_data::id)
       .def_readwrite("operating_code", &Fortran_mswell_data::operating_code)
       .def_readwrite("radius", &Fortran_mswell_data::radius)
       .def_readwrite("maximum_pressure",
                      &Fortran_mswell_data::maximum_pressure)
       .def_readwrite("minimum_pressure",
                      &Fortran_mswell_data::minimum_pressure)
       .def_readwrite("imposed_flowrate",
                      &Fortran_mswell_data::imposed_flowrate)
       .def_readwrite("injection_temperature",
                      &Fortran_mswell_data::injection_temperature)

       .def_readwrite("well_type", &Fortran_mswell_data::well_type)
       .def_property_readonly("is_closed",
                              [](const Fortran_mswell_data& self) {
                                 return self.operating_code == 'c';
                              })
       .def("open",
            [](Fortran_mswell_data& self) { self.operating_code = 'f'; })
       .def("close",
            [](Fortran_mswell_data& self) { self.operating_code = 'c'; });

   module.def("nb_mswells", &nb_mswells);

   module.def(
       "mswells_data",
       [](bool own_only) {
          auto p = get_mswells_data();
          const std::size_t n =
              own_only ? number_of_own_mswells() : nb_mswells();
          return py::make_iterator(p, p + n);
       },
       py::return_value_policy::reference, py::arg("own_only") = false);
}
