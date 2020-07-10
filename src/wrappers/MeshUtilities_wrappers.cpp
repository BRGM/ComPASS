//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

#include "ArrayWrapper.h"
#include "MeshUtilities.h"
#include "PyArrayWrapper.h"
#include "PyBuffer_wrappers.h"
#include "PyXArrayWrapper.h"
#include "XArrayWrapper.h"

// FUTURE: Code information bitwise?
struct NodeInfo {
   char proc;      // 'o' / 'g': own / ghost
   char frac;      // 'y' / 'n' : node in fracture / not in fracture
                   // FIXME: Name is to be changed
   char pressure;  // 'd' / 'n' / 'i' : dirichlet / neumann / interior for the
                   // pressure
                   // FIXME: preprocessor directives to be removed!
   char temperature;  // 'd' / 'n' / 'i' : dirichlet / neumann / interior for
                      // the temperature
};

// Fortran functions
#include "MeshSchema.fh"
extern "C" {
void get_global_number_of_nodes(long&);
void get_global_number_of_cells(long&);
void retrieve_global_vertices(XArrayWrapper<double>&);
void retrieve_vertices(XArrayWrapper<double>&);
void retrieve_cell_centers(XArrayWrapper<double>&);
void retrieve_face_centers(XArrayWrapper<double>&);
void retrieve_nodeflags(XArrayWrapper<int>&);
void retrieve_global_nodeflags(XArrayWrapper<int>&);
void retrieve_cellflags(XArrayWrapper<int>&);
void retrieve_global_cellflags(XArrayWrapper<int>&);
void retrieve_celltypes(XArrayWrapper<int8_t>&);
void retrieve_global_celltypes(XArrayWrapper<int8_t>&);
void retrieve_facetypes(XArrayWrapper<int8_t>&);
void retrieve_global_facetypes(XArrayWrapper<int8_t>&);
void retrieve_faceflags(XArrayWrapper<int>&);
void retrieve_global_faceflags(XArrayWrapper<int>&);
void retrieve_global_cell_rocktypes(XArrayWrapper<int>&);
void retrieve_global_node_rocktypes(XArrayWrapper<int>&);
void retrieve_global_fracture_rocktypes(XArrayWrapper<int>&);
void retrieve_cell_rocktypes(XArrayWrapper<int>&);
void retrieve_node_rocktypes(XArrayWrapper<int>&);
void retrieve_fracture_rocktypes(XArrayWrapper<int>&);
void retrieve_global_mesh_connectivity(MeshConnectivity&);
void retrieve_mesh_connectivity(MeshConnectivity&);
void retrieve_nodes_by_fractures(COC&);
void retrieve_global_id_faces(ArrayWrapper&);
void retrieve_global_cell_porosity(ArrayWrapper&);
void retrieve_global_fracture_porosity(ArrayWrapper&);
void retrieve_global_cell_permeability(ArrayWrapper&);
void retrieve_global_fracture_permeability(ArrayWrapper&);
void retrieve_cell_permeability(ArrayWrapper&);
void retrieve_fracture_permeability(ArrayWrapper&);
void retrieve_cell_porosity(ArrayWrapper&);
void retrieve_fracture_porosity(ArrayWrapper&);
#ifdef _THERMIQUE_
void retrieve_allthermalsources(XArrayWrapper<double>&);
void retrieve_global_cell_thermal_conductivity(ArrayWrapper&);
void retrieve_global_fracture_thermal_conductivity(ArrayWrapper&);
void retrieve_cell_thermal_conductivity(ArrayWrapper&);
void retrieve_fracture_thermal_conductivity(ArrayWrapper&);
void retrieve_cellthermalsource(XArrayWrapper<double>&);
void retrieve_nodethermalsource(XArrayWrapper<double>&);
void retrieve_fracthermalsource(XArrayWrapper<double>&);
void retrieve_cell_heat_source(ArrayWrapper&);
void retrieve_all_Fourier_porous_volumes(XArrayWrapper<double>&);
void retrieve_porovolfouriercell(XArrayWrapper<double>&);
void retrieve_porovolfouriernode(XArrayWrapper<double>&);
#endif
void retrieve_global_id_node(XArrayWrapper<NodeInfo>&);
void retrieve_id_node(XArrayWrapper<NodeInfo>&);
void retrieve_frac_face_id(XArrayWrapper<int>&);
void retrieve_face_frac_id(XArrayWrapper<int>&);
void retrieve_nb_cells_own(XArrayWrapper<int>&);
void retrieve_nb_faces_own(XArrayWrapper<int>&);
void retrieve_nb_nodes_own(XArrayWrapper<int>&);
void retrieve_nb_fractures_own(XArrayWrapper<int>&);
void retrieve_nb_wellinj_own(XArrayWrapper<int>&);
void retrieve_nb_wellprod_own(XArrayWrapper<int>&);
}

namespace py = pybind11;

template <typename Buffer, typename PyClass, typename Accessor>
void register_property_buffer(PyClass& pyclass, const char* name,
                              Accessor accessor) {
   pyclass.def_property_readonly(
       name,
       [accessor](py::object& self) {
          typedef typename Buffer::value_type value_type;
          auto buffer = retrieve_buffer<Buffer>(accessor).buffer_info();
          return py::array_t<value_type, py::array::c_style>{
              buffer.shape, buffer.strides,
              reinterpret_cast<value_type*>(buffer.ptr), self};
       },
       // py::keep_alive<0, 1>(): because we want the holder instance to be kept
       // alive as long as the buffer is used - should be ok this is the default
       // for def_property it's useless till the buffers are allocated globally
       // on the Fortran side but could become usefull later
       py::return_value_policy::reference_internal);
}

void add_mesh_utilities_wrappers(py::module& module) {
   add_vertices_array_wrapper(module, "global_vertices",
                              retrieve_global_vertices);
   add_vertices_array_wrapper(module, "vertices", retrieve_vertices);
   add_vertices_array_wrapper(module, "cell_centers", retrieve_cell_centers);
   add_vertices_array_wrapper(module, "face_centers", retrieve_face_centers);
   add_rocktypes_array_wrapper(module, "global_cell_rocktypes",
                               retrieve_global_cell_rocktypes);
   add_rocktypes_array_wrapper(module, "global_node_rocktypes",
                               retrieve_global_node_rocktypes);
   add_rocktypes_array_wrapper(module, "global_fracture_rocktypes",
                               retrieve_global_fracture_rocktypes);
   add_rocktypes_array_wrapper(module, "cell_rocktypes",
                               retrieve_cell_rocktypes);
   add_rocktypes_array_wrapper(module, "node_rocktypes",
                               retrieve_node_rocktypes);
   add_rocktypes_array_wrapper(module, "fracture_rocktypes",
                               retrieve_fracture_rocktypes);
   add_array_wrapper(module, "global_nodeflags", retrieve_global_nodeflags);
   add_array_wrapper(module, "global_cellflags", retrieve_global_cellflags);
   add_array_wrapper(module, "global_faceflags", retrieve_global_faceflags);
   add_array_wrapper(module, "global_celltypes", retrieve_global_celltypes);
   add_array_wrapper(module, "global_facetypes", retrieve_global_facetypes);
   add_array_wrapper(module, "nodeflags", retrieve_nodeflags);
   add_array_wrapper(module, "cellflags", retrieve_cellflags);
   add_array_wrapper(module, "faceflags", retrieve_faceflags);
   add_array_wrapper(module, "celltypes", retrieve_celltypes);
   add_array_wrapper(module, "facetypes", retrieve_facetypes);
   add_array_wrapper(module, "face_frac_id", retrieve_face_frac_id);
   add_array_wrapper(module, "frac_face_id", retrieve_frac_face_id);
   add_array_wrapper(module, "nb_cells_own", retrieve_nb_cells_own);
   add_array_wrapper(module, "nb_faces_own", retrieve_nb_faces_own);
   add_array_wrapper(module, "nb_nodes_own", retrieve_nb_nodes_own);
   add_array_wrapper(module, "nb_fractures_own", retrieve_nb_fractures_own);
   add_array_wrapper(module, "nb_wellinj_own", retrieve_nb_wellinj_own);
   add_array_wrapper(module, "nb_wellprod_own", retrieve_nb_wellprod_own);
#ifdef _THERMIQUE_
   add_array_wrapper(module, "all_thermal_sources", retrieve_allthermalsources);
   add_array_wrapper(module, "cellthermalsource", retrieve_cellthermalsource);
   add_array_wrapper(module, "nodethermalsource", retrieve_nodethermalsource);
   add_array_wrapper(module, "fracturethermalsource",
                     retrieve_fracthermalsource);
   add_array_wrapper(module, "all_Fourier_porous_volumes",
                     retrieve_all_Fourier_porous_volumes);
   add_array_wrapper(module, "porovolfouriercell", retrieve_porovolfouriercell);
   add_array_wrapper(module, "porovolfouriernode", retrieve_porovolfouriernode);
#endif

   module.def("global_number_of_nodes", []() {
      long n = 0;
      get_global_number_of_nodes(n);
      return n;
   });

   module.def("global_number_of_cells", []() {
      long n = 0;
      get_global_number_of_cells(n);
      return n;
   });

   module.def(
       "get_global_id_faces_buffer",
       []() { return retrieve_buffer<IntBuffer>(retrieve_global_id_faces); },
       "Get faces integer flag. Can be used to specify fracture faces setting "
       "the flag to -2.");

   module.def("get_global_cell_heat_source_buffer", []() {
      return retrieve_buffer<DoubleBuffer>(retrieve_cell_heat_source);
   });

   module.def("get_global_cell_porosity_buffer", []() {
      return retrieve_buffer<DoubleBuffer>(retrieve_global_cell_porosity);
   });

   module.def("get_global_fracture_porosity_buffer", []() {
      return retrieve_buffer<DoubleBuffer>(retrieve_global_fracture_porosity);
   });

   module.def("get_global_cell_permeability_buffer", []() {
      return retrieve_buffer<TensorBuffer>(retrieve_global_cell_permeability);
   });

   module.def("get_global_fracture_permeability_buffer", []() {
      return retrieve_buffer<DoubleBuffer>(
          retrieve_global_fracture_permeability);
   });

#ifdef _THERMIQUE_

   module.def("get_global_cell_thermal_conductivity_buffer", []() {
      return retrieve_buffer<TensorBuffer>(
          retrieve_global_cell_thermal_conductivity);
   });

   module.def("get_global_fracture_thermal_conductivity_buffer", []() {
      return retrieve_buffer<DoubleBuffer>(
          retrieve_global_fracture_thermal_conductivity);
   });

#endif

   // A holder class that could grow later
   // This is used so that numpy arrays are views and not copies
   struct Petrophysics {};

   auto petrophysics = py::class_<Petrophysics>(module, "petrophysics");
   petrophysics.def(py::init<>());
   register_property_buffer<DoubleBuffer>(petrophysics, "cell_porosity",
                                          retrieve_cell_porosity);
   register_property_buffer<TensorBuffer>(petrophysics, "cell_permeability",
                                          retrieve_cell_permeability);
   register_property_buffer<DoubleBuffer>(petrophysics, "fracture_porosity",
                                          retrieve_fracture_porosity);
   register_property_buffer<DoubleBuffer>(petrophysics, "fracture_permeability",
                                          retrieve_fracture_permeability);
#ifdef _THERMIQUE_
   register_property_buffer<TensorBuffer>(petrophysics,
                                          "cell_thermal_conductivity",
                                          retrieve_cell_thermal_conductivity);
   register_property_buffer<DoubleBuffer>(
       petrophysics, "fracture_thermal_conductivity",
       retrieve_fracture_thermal_conductivity);
#endif

   py::class_<MeshConnectivity>(module, "MeshConnectivity")
       .def_readwrite("NodebyCell", &MeshConnectivity::NodebyCell)
       .def_readwrite("NodebyFace", &MeshConnectivity::NodebyFace)
       .def_readwrite("FacebyNode", &MeshConnectivity::FacebyNode)
       .def_readwrite("FacebyCell", &MeshConnectivity::FacebyCell)
       .def_readwrite("CellbyNode", &MeshConnectivity::CellbyNode)
       .def_readwrite("CellbyFace", &MeshConnectivity::CellbyFace)
       .def_readwrite("CellbyCell", &MeshConnectivity::CellbyCell);

   module.def(
       "get_global_connectivity",
       []() {
          MeshConnectivity connectivity;
          retrieve_global_mesh_connectivity(connectivity);
          return connectivity;
       },
       "Get global mesh connectivity.");

   module.def(
       "get_connectivity",
       []() {
          MeshConnectivity connectivity;
          retrieve_mesh_connectivity(connectivity);
          return connectivity;
       },
       "Get local mesh connectivity.");

   module.def("get_nodes_by_fractures", []() {
      COC coc;
      retrieve_nodes_by_fractures(coc);
      return coc;
   });

   PYBIND11_NUMPY_DTYPE(NodeInfo, proc, frac, pressure, temperature);

   add_array_wrapper(module, "global_node_info", retrieve_global_id_node);
   add_array_wrapper(module, "node_info", retrieve_id_node);

   module.def("number_of_nodes", &number_of_nodes);
   module.def("number_of_own_nodes", &number_of_own_nodes);
   module.def("number_of_cells", &number_of_cells);
   module.def("number_of_own_cells", &number_of_own_cells);
   module.def("number_of_fractures", &number_of_fractures);
   module.def("number_of_own_fractures", &number_of_own_fractures);
   module.def("number_of_own_injectors", &number_of_own_injectors);
   module.def("number_of_own_producers", &number_of_own_producers);
}
