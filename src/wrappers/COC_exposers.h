//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <limits>
#include <string>

#include "COC.h"

namespace py = pybind11;

template <typename Value_type, bool with_buffer_protocol>
struct COC_container_factory;

template <typename Value_type>
struct COC_container_factory<Value_type, true> {
   static auto make(py::module& module, std::string basename) {
      return py::class_<GenericCOC_container<Value_type>>(
          module, basename.append("COCcontainer").c_str(),
          py::buffer_protocol());
   }
};

template <typename Value_type>
struct COC_container_factory<Value_type, false> {
   static auto make(py::module& module, std::string basename) {
      return py::class_<GenericCOC_container<Value_type>>(
          module, basename.append("COCcontainer").c_str());
   }
};

template <typename Value_type, bool with_buffer_protocol = false>
auto expose_coc(py::module& module, std::string basename) {
   typedef GenericCOC_container<Value_type> Coc_cont;
   typedef GenericCOC_iterator<Value_type> Coc_it;
   typedef GenericCOC<Value_type> Coc;

   auto container =
       COC_container_factory<Value_type, with_buffer_protocol>::make(module,
                                                                     basename);
   container.def("__len__",
                 [](const COC_container& cocc) { return cocc.length(); });

   auto iterator =
       py::class_<Coc_it>(module, basename.append("COCiterator").c_str());

   auto coc =
       py::class_<Coc>(module, basename.append("COC").c_str())
           .def(
               "__iter__",
               [](Coc& coc) {
                  return py::make_iterator(coc.begin(), coc.end());
               },
               py::keep_alive<0,
                              1>() /* Keep Coc alive while iterator is used */
               )
           .def("__getitem__",
                [](Coc& coc, int i) -> Coc_cont { return coc[i]; })
           .def("__len__",
                [](const Coc& coc) { return coc.number_of_containers(); })
           .def(
               "offsets",
               [](py::object object) {
                  auto coc = object.cast<Coc&>();
                  assert(coc.number_of_containers() + 1 <
                         std::numeric_limits<py::ssize_t>::max());
                  return py::array_t<typename Coc::offset_type,
                                     py::array::c_style>{
                      {static_cast<py::ssize_t>(coc.number_of_containers() +
                                                1)},
                      coc.offset_data(),
                      object};
               },
               py::keep_alive<0, 1>())
           .def(
               "contiguous_content",
               [](py::object object) {
                  auto coc = object.cast<Coc&>();
                  assert(coc.size() < std::numeric_limits<py::ssize_t>::max());
                  return py::array_t<typename Coc::value_type,
                                     py::array::c_style>{
                      {static_cast<py::ssize_t>(coc.size())},
                      coc.content_data(),
                      object};
               },
               py::keep_alive<0, 1>());

   return std::make_tuple(coc, iterator, container);
}
