//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/pybind11.h>

#include "COC.h"

struct DOFId {
   int proc;
   int local_id;
};

typedef GenericCOC<DOFId> DofFamilyCOC;

// Fortran functions
extern "C" {
void retrieve_NumNodebyProc(DofFamilyCOC&);
void retrieve_NumFracbyProc(DofFamilyCOC&);
void retrieve_NumWellProdbyProc(DofFamilyCOC&);
void retrieve_NumWellInjbyProc(DofFamilyCOC&);
}

#include "COC_exposers.h"

namespace py = pybind11;

void add_mesh_schema_wrappers(py::module& module) {
   typedef GenericCOC_container<DOFId> DofFamily;

   auto exposed = expose_coc<DOFId>(module, "DofFamily_");
   auto& coc_class = std::get<0>(exposed);
   auto& container_class = std::get<2>(exposed);

   container_class.def("as_array", [](DofFamily& self) {
      std::array<std::size_t, 2> shape{
          {self.length(), static_cast<std::size_t>(2)}};
      std::array<std::size_t, 2> stride{{2 * sizeof(int), sizeof(int)}};
      return py::array_t<int, py::array::c_style>{
          shape, stride, reinterpret_cast<int*>(self.data())};
   });
   container_class.def_property_readonly("proc", [](DofFamily& self) {
      std::array<std::size_t, 1> shape{{self.length()}};
      std::array<std::size_t, 1> stride{{2 * sizeof(int)}};
      return py::array_t<int, py::array::c_style>{
          shape, stride, reinterpret_cast<int*>(self.data())};
   });
   container_class.def_property_readonly("local_id", [](DofFamily& self) {
      std::array<std::size_t, 1> shape{{self.length()}};
      std::array<std::size_t, 1> stride{{2 * sizeof(int)}};
      return py::array_t<int, py::array::c_style>{
          shape, stride, reinterpret_cast<int*>(self.data()) + 1};
   });

   coc_class.def("number_of_domains", &DofFamilyCOC::number_of_containers);
   coc_class.def("as_array", [](DofFamilyCOC& self) {
      std::array<std::size_t, 2> shape{
          {self.size(), static_cast<std::size_t>(2)}};
      std::array<std::size_t, 2> stride{{2 * sizeof(int), sizeof(int)}};
      return py::array_t<int, py::array::c_style>{
          shape, stride, reinterpret_cast<int*>(self.mutable_content_data())};
   });

   auto expose = [&module](const char* name, void (*retrieve)(DofFamilyCOC&)) {
      module.def(name, [retrieve]() {
         auto coc = std::make_unique<DofFamilyCOC>();
         (*retrieve)(*coc);
         return coc;
      });
   };

   expose("NumNodebyProc", retrieve_NumNodebyProc);
   expose("NumFracbyProc", retrieve_NumFracbyProc);
   expose("NumWellInjbyProc", retrieve_NumWellInjbyProc);
   expose("NumWellProdbyProc", retrieve_NumWellProdbyProc);
}
