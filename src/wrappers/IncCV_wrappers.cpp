//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <array>
#include <cstddef>
#include <sstream>

#include "Neumann.h"
#include "StateObjects.h"

// Fortran functions
extern "C" {
int number_of_components();
int number_of_phases();
int size_of_unknowns();
void retrieve_dirichlet_node_states(StateArray&);
void retrieve_all_states(StateArray&);
void retrieve_node_states(StateArray&);
void retrieve_fracture_states(StateArray&);
void retrieve_cell_states(StateArray&);
void retrieve_own_dirichlet_node_states(StateArray&);
void retrieve_own_node_states(StateArray&);
void retrieve_own_fracture_states(StateArray&);
void retrieve_own_cell_states(StateArray&);
void retrieve_injection_whp(XArrayWrapper<double>&);
void retrieve_production_whp(XArrayWrapper<double>&);
void retrieve_production_whp(XArrayWrapper<double>&);
void set_face_with_neumann_contribution(int, const NeumannBC&);
void set_faces_with_neumann_contribution(int, const int*, const NeumannBC&);
void set_face_neumann_contributions(int, int*, NeumannBC*);
void set_fracture_edge_with_neumann_contribution(const int*, const NeumannBC&);
void set_fracture_edges_with_neumann_contribution(int, const int*,
                                                  const NeumannBC&);
void set_fracture_edge_neumann_contributions(int, int*, NeumannBC*);
void clear_all_neumann_contributions();
void LoisThermoHydro_total_specific_enthalpy(const X&, double& h);
#ifndef NDEBUG
void dump_incv_info();
#endif
}

#include <pybind11/numpy.h>

#include "PyXArrayWrapper.h"
#include "attribute_array.h"
#include "wrap_array_holder.h"

namespace py = pybind11;

void add_IncCV_wrappers(py::module& module) {
   using Context = typename X::Model::Context;
   using Real = typename X::Model::Real;
   constexpr auto np = X::Model::np;  // number of phases // FIXME to be passed
                                      // dynamically as argument
   constexpr auto nc = X::Model::nc;  // number of components // FIXME to be
                                      // passed dynamically as argument
   constexpr auto nbdof = X::Model::nbdof;  // number of components // FIXME to
                                            // be passed dynamically as argument

   module.def("check_IncCV", []() {
      py::print("Check IncCV:", number_of_components(), number_of_phases(),
                size_of_unknowns(), sizeof(X));
   });

   py::class_<X>(module, "State")
       .def(py::init<const X&>())
       .def_readwrite("context", &X::context)
       .def_readwrite("p", &X::p)
       .def_property_readonly("pa",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     np, self.cast<X&>().pa.data(), self};
                              })
       .def_readwrite("T", &X::T)
       .def_property_readonly("S",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     np, self.cast<X&>().S.data(), self};
                              })
       .def_property_readonly(
           "C",
           [](py::object& self) {
              return py::array_t<Real, py::array::c_style>{
                  {np, nc}, self.cast<X&>().C.data()->data(), self};
           })
       .def_property_readonly("accumulation",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     nc, self.cast<X&>().accumulation.data(),
                                     self};
                              })
       .def("__str__", [](const X& self) {
          std::basic_stringstream<char> s;
          s << "Context: " << static_cast<int>(self.context);
          s << " p=" << self.p;
          s << " T=" << self.T;
          s << " S=(";
          for (auto&& Sk : self.S) s << Sk << ", ";
          s << ") C=(";
          for (auto&& Ck : self.C) {
             s << "(";
             for (auto&& Cki : Ck) s << Cki << ", ";
             s << "), ";
          }
          s << ")";
          return py::str{s.str()};
       });

   module.def("build_state",
              [](const X& other) { return std::make_unique<X>(other); });

   // FUTURE: all the following until next FUTURE tag shall be useless soon (cf
   // infra)
   auto pyStateArray = wrap_array_holder<StateArray>(module, "States");
   add_attribute_array<StateArray, Context>(pyStateArray, "context",
                                            offsetof(X, context));
   add_attribute_array<StateArray, Real>(pyStateArray, "p", offsetof(X, p));
   add_attribute_array<StateArray, Real>(pyStateArray, "pa", offsetof(X, pa),
                                         {np}, {sizeof(Real)});
   add_attribute_array<StateArray, Real>(pyStateArray, "T", offsetof(X, T));
   add_attribute_array<StateArray, Real>(pyStateArray, "C", offsetof(X, C),
                                         {np, nc},
                                         {nc * sizeof(Real), sizeof(Real)});
   add_attribute_array<StateArray, Real>(pyStateArray, "S", offsetof(X, S),
                                         {np}, {sizeof(Real)});
   add_attribute_array<StateArray, Real>(pyStateArray, "accumulation",
                                         offsetof(X, accumulation), {nbdof},
                                         {sizeof(Real)});
#ifdef _WIP_FREEFLOW_STRUCTURES_
   add_attribute_array<StateArray, Real>(
       pyStateArray, "FreeFlow_phase_flowrate",
       offsetof(X, FreeFlow_phase_flowrate), {np}, {sizeof(Real)});
#endif  // _WIP_FREEFLOW_STRUCTURES_
   // auto PyNeumannArray = py::class_<NeumannArray>(module,
   // "NeumannContributions")
   //    .def("size", [](const NeumannArray& self) { return self.length; })
   //    .def_property_readonly("shape", [](const NeumannArray& self) { return
   //    py::make_tuple(self.length); });
   // add_attribute_array<typename Neumann::molar_flux>(PyNeumannArray,
   // "molar_flux", offsetof(Neumann, molar_flux), { Neumann::nc }, {
   // sizeof(Real) }); add_attribute_array<typename
   // Neumann::Real>(PyNeumannArray, "heat_flux", offsetof(Neumann, heat_flux));

   module.def("dirichlet_node_states", []() {
      return StateArray::retrieve(retrieve_dirichlet_node_states);
   });
   module.def(
       "all_states", []() { return StateArray::retrieve(retrieve_all_states); },
       R"doc(
        Retrieve all nodes, fractures and cell states in a contiguous structure.
    )doc");
   module.def("node_states",
              []() { return StateArray::retrieve(retrieve_node_states); });
   module.def("fracture_states",
              []() { return StateArray::retrieve(retrieve_fracture_states); });
   module.def("cell_states",
              []() { return StateArray::retrieve(retrieve_cell_states); });
   module.def("own_dirichlet_node_states", []() {
      return StateArray::retrieve(retrieve_own_dirichlet_node_states);
   });
   module.def("own_node_states",
              []() { return StateArray::retrieve(retrieve_own_node_states); });
   module.def("own_fracture_states", []() {
      return StateArray::retrieve(retrieve_own_fracture_states);
   });
   module.def("own_cell_states",
              []() { return StateArray::retrieve(retrieve_own_cell_states); });

   module.def("number_of_components",
              []() { return ComPASS_NUMBER_OF_COMPONENTS; });
   module.def("number_of_phases", []() { return ComPASS_NUMBER_OF_PHASES; });

   // FUTURE: PYBIND11_NUMPY_DTYPE shall support arrays in forthcoming release
   // of pybind11 PYBIND11_NUMPY_DTYPE(X, context, p, T, C, S, accumulation);
   // Yet with a slightly cumbersome syntax: X['p'] to access pressure instead
   // of X.p
   // PYBIND11_NUMPY_DTYPE(X, context, p, T);
   // module.def("node_states", []() {
   //	auto wrapper = StateArray{};
   //	retrieve_boundary_states(wrapper);
   //	return py::array_t<X, py::array::c_style>{wrapper.length,
   // wrapper.pointer};
   //});

   static_assert(
       std::is_same<typename X::Model, typename NeumannBC::Model>::value,
       "Model inconsistency");

   py::class_<NeumannBC>(module, "NeumannBC")
       .def(py::init<>())
       .def_readwrite("heat_flux", &NeumannBC::heat_flux)
       .def_property_readonly(
           "molar_flux",
           [](py::object& self) {
              auto data =
                  self.cast<NeumannBC*>()
                      ->molar_flux
                      .data();  // C++ pointer to underlying NeumannBC instance
              return py::array_t<Real, py::array::c_style>({nc}, {sizeof(Real)},
                                                           data, self);
           })  // return value policy defaults to
               // py::return_value_policy::reference_internal which is ok
       ;

   module.def(
       "set_Neumann_faces",
       [](py::array_t<int, py::array::c_style | py::array::forcecast> faces,
          const NeumannBC& condition) {
          if (faces.ndim() != 1 || faces.size() == 0) return;
          auto raw_faces = faces.unchecked<1>();
          set_faces_with_neumann_contribution(raw_faces.shape(0),
                                              raw_faces.data(0), condition);
       });

   module.def(
       "set_Neumann_fracture_edges",
       [](py::array_t<int, py::array::c_style | py::array::forcecast> edges,
          const NeumannBC& condition) {
          if (edges.ndim() != 2 || edges.size() == 0) return;
          auto raw_edges = edges.unchecked<2>();
          set_fracture_edges_with_neumann_contribution(
              raw_edges.shape(0), raw_edges.data(0, 0), condition);
       });

   module.def("clear_all_neumann_contributions",
              &clear_all_neumann_contributions);

   // PYBIND11_NUMPY_DTYPE(NeumannBC, molar_flux, heat_flux);
   // module.def("neumann_conditions", [](std::size_t n) {
   //    return py::array_t<NeumannBC, py::array::c_style>{ n };
   //});

   // module.def("node_states", []() {
   //	auto wrapper = StateArray{};
   //	retrieve_boundary_states(wrapper);
   //	return py::array_t<X, py::array::c_style>{wrapper.length,
   // wrapper.pointer};
   //});

   add_array_wrapper(module, "injection_whp", retrieve_injection_whp);
   add_array_wrapper(module, "production_whp", retrieve_production_whp);

#ifndef NDEBUG
   module.def("dump_incv_info", &dump_incv_info);
#endif

   module.def("total_specific_enthalpy", [](const X& x) {
      double h = 0;
      LoisThermoHydro_total_specific_enthalpy(x, h);
      return h;
   });
   module.def(
       "total_specific_enthalpy",
       [](const StateArray& states) {
          auto res = py::array_t<double, py::array::c_style>{
              static_cast<py::ssize_t>(states.length)};
          auto* p = res.mutable_data();
          for (auto&& X : states) {
             LoisThermoHydro_total_specific_enthalpy(X, *p);
             ++p;
          }
          return res;
       },
       py::return_value_policy::take_ownership);
}
