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

namespace py = pybind11;

#include "XArrayWrapper.h"

template <typename T>
inline auto encapsulate(T* p) {
   return py::capsule(p);
}

// FIXME: Creating local reference through py::array_t<typename
// Wrapper::wrapped_type, py::array::c_style>{} is a bit weird/dangerous ? User
// buffer_info instead or a DummyStructure so that the returned array is not a
// copy ?
template <typename Wrapper>
static auto retrieve_ndarray(std::function<void(Wrapper&)> bind) {
   auto wrapper = Wrapper{};
   bind(wrapper);
   assert(wrapper.length < std::numeric_limits<py::ssize_t>::max());
   return py::array_t<typename Wrapper::wrapped_type, py::array::c_style>{
       static_cast<py::ssize_t>(wrapper.length), wrapper.pointer,
       py::array_t<typename Wrapper::wrapped_type, py::array::c_style>{}};
}

// FUTURE: This is to be removed with C++17 and template deductions
template <typename Wrapper>
static auto retrieve_ndarray(void (*bind)(Wrapper&)) {
   return retrieve_ndarray(std::function<void(Wrapper&)>{bind});
}

template <typename Wrapper>
py::module& add_array_wrapper(py::module& module, const char* getter_name,
                              void (*bind)(Wrapper&)) {
   module.def(getter_name, [bind]() { return retrieve_ndarray(bind); });
   return module;
}

// FIXME: Creating local reference through py::array_t<typename
// Wrapper::wrapped_type, py::array::c_style>{} is a bit weird/dangerous ? User
// buffer_info instead or a DummyStructure so that the returned array is not a
// copy ?
template <typename Wrapper>
py::module& add_vertices_array_wrapper(py::module& module,
                                       const char* getter_name,
                                       void (*bind)(Wrapper&)) {
   static_assert(std::is_same<double, typename Wrapper::wrapped_type>::value,
                 "coordinates should be double");
   module.def(getter_name, [bind]() {
      auto wrapper = Wrapper{};
      bind(wrapper);
      assert(wrapper.length < std::numeric_limits<py::ssize_t>::max());
      return py::array_t<double, py::array::c_style>{
          {static_cast<py::ssize_t>(wrapper.length),
           static_cast<py::ssize_t>(3)},
          wrapper.pointer,
          encapsulate(wrapper.pointer)};
   });
   return module;
}

// FIXME: Creating local reference through py::array_t<typename
// Wrapper::wrapped_type, py::array::c_style>{} is a bit weird/dangerous ? User
// buffer_info instead or a DummyStructure so that the returned array is not a
// copy ?
template <typename Wrapper>
py::module& add_rocktypes_array_wrapper(py::module& module,
                                        const char* getter_name,
                                        void (*bind)(Wrapper&)) {
   static_assert(std::is_same<int, typename Wrapper::wrapped_type>::value,
                 "rocktypes should be integers");
   module.def(getter_name, [bind]() {
      auto wrapper = Wrapper{};
      bind(wrapper);
      assert(wrapper.length < std::numeric_limits<py::ssize_t>::max());
      return py::array_t<int, py::array::c_style>{
#ifdef _THERMIQUE_
          {static_cast<py::ssize_t>(wrapper.length),
           static_cast<py::ssize_t>(2)},
#else
            {static_cast<py::ssize_t>(wrapper.length)},
#endif
          wrapper.pointer,
          encapsulate(wrapper.pointer)};
   });
   return module;
}
