#include <pybind11/pybind11.h>

#include <iostream>
#include <string>

struct StringWrapper {
   const char* pointer;
   const size_t length;
   StringWrapper() = delete;
   StringWrapper(const StringWrapper&) = delete;
   void operator=(const StringWrapper&) = delete;
   StringWrapper(const std::string& s)
       : pointer(s.c_str()), length(s.length()) {}
};

void cpp_dump(const std::string& s) {
   std::cout << "C++ writes: " << s << std::endl;
}

extern "C" {
void fortran_dump(const StringWrapper&);
}

namespace py = pybind11;

PYBIND11_PLUGIN(passingstring) {
   py::module m("passingstring",
                "pybind11 example plugin - passing string to Fortran");

   m.def("cpp_dump", &cpp_dump, "A function which prints a string from C++.");
   m.def(
       "fortran_dump", [](const std::string& s) { fortran_dump(s); },
       "A function which prints a string from Fortran.");

   return m.ptr();
}
