#include "PyBuffer_wrappers.h"

void add_pybuffer_wrappers(py::module& module)
{

	IntBuffer::add_buffer_class(module, "IntBuffer");
	DoubleBuffer::add_buffer_class(module, "DoubleBuffer");
	CoordinatesBuffer::add_buffer_class(module, "CoordinatesBuffer");
	TensorBuffer::add_buffer_class(module, "TensorBuffer");

}