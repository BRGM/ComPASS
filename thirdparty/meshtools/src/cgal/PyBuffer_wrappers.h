#pragma once

#include <functional>
#include <vector>
#include <pybind11/pybind11.h>

#include "ArrayWrapper.h"

namespace py = pybind11;

template <typename T, std::size_t ...extend>
class PyBufferWrapper
{
protected:
	ArrayWrapper wrapper;
	std::size_t ndim;
	std::vector<std::size_t> shape;
	std::vector<std::size_t> strides;
	PyBufferWrapper() = delete;
public:
	explicit PyBufferWrapper(const ArrayWrapper& AW) :
		wrapper{ AW },
		ndim{ sizeof...(extend)+1 },
		shape{ wrapper.length, extend... },
		strides{} {
		assert(ndim > 0);
		assert(shape.size() == ndim);
		strides.resize(ndim);
		strides[ndim - 1] = sizeof(T);
		for (std::size_t k = 0; k < ndim - 1; ++k) {
			strides[ndim - 2 - k] = shape[ndim - 1 - k] * strides[ndim - 1 - k];
		}
	}
	auto buffer_info() -> py::buffer_info {
		return py::buffer_info(
			static_cast<T *>(wrapper.pointer),
			sizeof(T), py::format_descriptor<T>::format(),
			ndim, shape, strides
		);
	}
	static auto add_buffer_class(py::module& module, const std::string& classname) {
		typedef PyBufferWrapper<T, extend...> Buffer;
		return py::class_<Buffer>(module, classname.c_str(), py::buffer_protocol())
			.def_buffer([](Buffer& buffer) -> py::buffer_info { return buffer.buffer_info();  });
	}
};
