#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include "XArrayWrapper.h"


// FIXME: Creating local reference through py::array_t<typename Wrapper::wrapped_type, py::array::c_style>{}
// is a bit weird/dangerous ?
// User buffer_info instead or a DummyStructure so that the returned array is not a copy ?
template <typename Wrapper>
static auto retrieve_ndarray(std::function <void(Wrapper&)> bind)
{
	auto wrapper = Wrapper{};
	bind(wrapper);
	return py::array_t<typename Wrapper::wrapped_type, py::array::c_style>{
		wrapper.length, wrapper.pointer,
			py::array_t<typename Wrapper::wrapped_type, py::array::c_style>{} 
	};
}

// FUTURE: This is to be removed with C++17 and template deductions
template <typename Wrapper>
static auto retrieve_ndarray(void (*bind)(Wrapper&))
{
	return retrieve_ndarray(std::function <void(Wrapper&)>{ bind });
}

template <typename Wrapper>
py::module& add_array_wrapper(py::module& module, const char *getter_name, void(*bind)(Wrapper&))
{
	module.def(getter_name, [bind]() { return retrieve_ndarray(bind); });
	return module;
}