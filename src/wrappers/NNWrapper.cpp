#include <cstddef>

#include <StringWrapper.h>

#include <pybind11/pybind11.h>

struct Vertices
{
  double * p;
  std::size_t nb_points;
  Vertices():
    p(NULL),
    nb_points(0)
  {}
};

extern "C"
{
  void NN_init(const StringWrapper&, const StringWrapper&, const StringWrapper&);
  void NN_init_up_to_mesh(const StringWrapper&, const StringWrapper&, const StringWrapper&);
  void NN_init_phase2(const StringWrapper&);
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
  module.def("init_up_to_mesh", 
    [](const std::string& MeshFile, const std::string& LogFile, const std::string& OutputDir) { 
      NN_init_up_to_mesh(MeshFile, LogFile, OutputDir); 
    },
    "Initialisation of ComPASS up to the point where the global mesh is loaded.");
  module.def("init_phase2", 
    [](const std::string& OutputDir) { 
      NN_init_phase2(OutputDir); 
    },
    "Initialisation of ComPASS phase 2 : partition and distribute.");
  module.def("main_loop", [](int TimeIter, const std::string& OutputDir) { NN_main(TimeIter, OutputDir); },
    "Main loop of ComPASS.");
  module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");

  py::class_<Vertices>(module, "VerticesHandle", py::buffer_protocol())
     .def_buffer([](Vertices &V) -> py::buffer_info {
        return py::buffer_info(
          V.p,                               /* Pointer to buffer */
          sizeof(double),                          /* Size of one scalar */
          py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
          2,                                      /* Number of dimensions */
          { 3, V.nb_points },                 /* Buffer dimensions */
          { sizeof(double),             /* Strides (in bytes) for each index */
            sizeof(double) * 3 }
        );
      });
  module.def("get_vertices", []() {
      Vertices V;
      vertices_buffer(V);
      return V;
    },
    "Get node coordinates."
  );

  return module.ptr();

}
