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

#include <array>
#include <cstddef>
#include <sstream>

#include "ArrayWrapper.h"
#include "PyXArrayWrapper.h"
#include "StateObjects.h"
#include "XArrayWrapper.h"
#include "attribute_array.h"
#include "wrap_array_holder.h"

extern "C" {
void reset_freeflow_nodes();
}
#ifdef _WITH_FREEFLOW_STRUCTURES_
// Fortran functions
extern "C" {
void clear_freeflow_faces();
void set_atm_temperature(double);
void set_rain_temperature(double);
void set_atm_pressure(double);
void set_atm_rain_flux(double);
void set_freeflow_faces(const ArrayWrapper&);
void retrieve_freeflow_nodes_mask(XArrayWrapper<bool>&);
void retrieve_freeflow_nodes_area(XArrayWrapper<double>&);
void retrieve_freeflow_node_states(StateFFArray&);
}
#endif
namespace py = pybind11;

void add_freeflow_wrappers(py::module& module) {
   using Real = typename XFF::Model::Real;
   constexpr auto np = XFF::Model::np;  // number of phases // FIXME to be
                                        // passed dynamically as argument
   constexpr auto nc = XFF::Model::nc;  // number of components // FIXME to be
                                        // passed dynamically as argument

   py::class_<XFF>(module, "FarFieldState")
       .def(py::init<const XFF&>())
       .def_readwrite("p", &XFF::p)
       .def_property_readonly("T",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     np, self.cast<XFF&>().T.data(), self};
                              })
       .def_property_readonly(
           "C",
           [](py::object& self) {
              return py::array_t<Real, py::array::c_style>{
                  {np, nc}, self.cast<XFF&>().C.data()->data(), self};
           })
       .def_property_readonly("imposed_flux",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     np, self.cast<XFF&>().imposed_flux.data(),
                                     self};
                              })
       .def_property_readonly("Hm",
                              [](py::object& self) {
                                 return py::array_t<Real, py::array::c_style>{
                                     np, self.cast<XFF&>().Hm.data(), self};
                              })
       .def_readwrite("HT", &XFF::HT)
       .def("__str__", [](const XFF& self) {
          std::basic_stringstream<char> s;
          s << " p=" << self.p;
          s << " T=(";
          for (auto&& Tk : self.T) s << Tk << ", ";
          s << ") C=(";
          for (auto&& Ck : self.C) {
             s << "(";
             for (auto&& Cki : Ck) s << Cki << ", ";
             s << "), ";
          }
          s << ") imposed_flux=(";
          for (auto&& Fk : self.imposed_flux) s << Fk << ", ";
          s << ") Hm (all phases)=(";
          for (auto&& Hmk : self.Hm) s << Hmk << ", ";
          s << ") HT=" << self.HT;
          return py::str{s.str()};
       });

   // FUTURE: all the following until next FUTURE tag shall be useless soon (cf
   // infra)
   auto pyStateFFArray =
       wrap_array_holder<StateFFArray>(module, "FarFieldStates");
   add_attribute_array<XFF, StateFFArray, Real>(pyStateFFArray, "p",
                                                offsetof(XFF, p));
   add_attribute_array<XFF, StateFFArray, Real>(
       pyStateFFArray, "T", offsetof(XFF, T), {np}, {sizeof(Real)});
   add_attribute_array<XFF, StateFFArray, Real>(
       pyStateFFArray, "C", offsetof(XFF, C), {np, nc},
       {nc * sizeof(Real), sizeof(Real)});
   add_attribute_array<XFF, StateFFArray, Real>(pyStateFFArray, "imposed_flux",
                                                offsetof(XFF, imposed_flux),
                                                {np}, {sizeof(Real)});
   add_attribute_array<XFF, StateFFArray, Real>(
       pyStateFFArray, "Hm", offsetof(XFF, Hm), {np}, {sizeof(Real)});
   add_attribute_array<XFF, StateFFArray, Real>(pyStateFFArray, "HT",
                                                offsetof(XFF, HT));

#ifdef _WITH_FREEFLOW_STRUCTURES_
   module.attr("has_freeflow_structures") = true;
#else
   module.attr("has_freeflow_structures") = false;
#endif

#ifdef _WITH_FREEFLOW_STRUCTURES_
   module.def("clear_freeflow_faces", &clear_freeflow_faces);
   module.def("reset_freeflow_nodes", &reset_freeflow_nodes);
   module.def("set_atm_temperature", &set_atm_temperature);
   module.def("set_rain_temperature", &set_rain_temperature);
   module.def("set_atm_pressure", &set_atm_pressure);
   module.def("set_atm_rain_flux", &set_atm_rain_flux);
   module.def(
       "set_freeflow_faces",
       [](py::array_t<int, py::array::c_style | py::array::forcecast> faces) {
          if (faces.ndim() != 1)
             throw std::runtime_error(
                 "Freeflow faces array shoud be a single dimension array.");
          if (faces.size() == 0) return;  // empty list, does nothing
          auto wrapper =
              ArrayWrapper::wrap(faces.mutable_data(0), faces.size());
          set_freeflow_faces(wrapper);
       });
   module.def("retrieve_freeflow_nodes_mask",
              []() { return retrieve_ndarray(retrieve_freeflow_nodes_mask); });
   module.def("retrieve_freeflow_nodes_area",
              []() { return retrieve_ndarray(retrieve_freeflow_nodes_area); });
   module.def("freeflow_node_states", []() {
      return StateFFArray::retrieve(retrieve_freeflow_node_states);
   });
#else
   module.def("reset_freeflow_nodes", []() {});
#endif
}
