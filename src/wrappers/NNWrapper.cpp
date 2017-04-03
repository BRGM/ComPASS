#include <cstddef>

#include <StringWrapper.h>

#include <pybind11/pybind11.h>

struct Vertices
{
  double * p;
  std::size_t n;
};

extern "C"
{
  void NN_init(const StringWrapper&, const StringWrapper&, const StringWrapper&);
  void NN_main(int, const StringWrapper&);
  void NN_finalize();
  void vertices_buffer(Vertices&);
}

namespace py = pybind11;

PYBIND11_PLUGIN(ComPASS)
{

  py::module module("ComPASS", "pybind11 ComPASS library interface");

  module.def("init", 
    [](const std::string& MeshFile, const std::string& LogFile, const std::string& OutputDir) { 
      NN_init(MeshFile, LogFile, OutputDir); 
    },
    "Initialisation of ComPASS.");
  module.def("main_loop", [](int TimeIter, const std::string& OutputDir) { NN_main(TimeIter, OutputDir); },
    "Main loop of ComPASS.");
  module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");

  module.def("get_vertices", []() {Vertices V; vertices_buffer(V); }, "Get node coordinates.");

  return module.ptr();

}
