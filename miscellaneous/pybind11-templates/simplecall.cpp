#include <pybind11/pybind11.h>

int add(int i, int j)
{
	return i + j;
}


extern "C"
{
	int my_fortran_add(int, int);
	int my_fortran_add_byref(int&, int&);
}

namespace py = pybind11;

PYBIND11_PLUGIN(simplecall)
{

	py::module m("simplecall", "pybind11 example plugin - simple function call");

	m.def("add", &add, "A function which adds two numbers in C++.");
	m.def("fadd", &my_fortran_add, "A function which adds two numbers in Fortran.");
	m.def("fadd_byref", &my_fortran_add_byref, "A function which adds two numbers in Fortran.");

	return m.ptr();

}
