//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "SyncPetsc_wrappers.h"

#include <petscmat.h>
#include <petscvec.h>

#include "NewtonIncrements.h"

// Passing PETSc objects between C and Fortran does not rely on iso C binding
// but on opaque pointers and compiler dependant name mangling hence some
// restrictions on the case of the subroutine names. It is supposed that the
// compiler expose the Fortran routine with lower case and a trailing underscore
// (what gcc does). In the above PETSc example a more elaborated strategy is
// used (based on multiple preprocessor defines... !). As this wrapping is bound
// to disappear (one day) we keep it "simple" but not very portable.

// Fortran functions
extern "C" {
void syncpetsc_getsolnodefracwellmswell_(Vec*,
                                         NewtonIncrements::Pointers<double>);
void syncpetsc_getsolnodefracwell_(Vec*, NewtonIncrements::Pointers<double>);
void SyncPetsc_colnum(int*, std::size_t);
void MeshSchema_part_info(PartInfo&);
void SyncPetscMSWells_colnum(int*, std::size_t);
void syncpetscmswells_getsolmswell_(Vec*, NewtonIncrements::Pointers<double>);
}

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Petsc_caster.h"

template <typename T>
auto retrieve_matrix_size(T accessor) {
   auto res = py::array_t<MatrixSize::size_type, py::array::c_style>{2};
   accessor(*reinterpret_cast<MatrixSize*>(res.mutable_data()));
   return res;
}

void add_SyncPetsc_wrappers(py::module& module) {
   module.def("SyncPetsc_GetSolNodeFracWell",
              [](py::object V, NewtonIncrements& increments) {
                 auto vec = cast_to_PETSc<Vec>(V);
                 syncpetsc_getsolnodefracwell_(&vec, increments.pointers());
              });

   module.def("SyncPetsc_GetSolNodeFracWellMSWell",
              [](py::object V, NewtonIncrements& increments) {
                 auto vec = cast_to_PETSc<Vec>(V);
                 syncpetsc_getsolnodefracwellmswell_(&vec,
                                                     increments.pointers());
              });

   //     module.def("SyncPetsc_global_matrix_size", []() {
   //           return retrieve_matrix_size(SyncPetsc_global_matrix_size);
   //     });

   //     module.def("SyncPetsc_local_matrix_size", []() {
   //           return retrieve_matrix_size(SyncPetsc_local_matrix_size);
   //     });

   py::class_<PartElement>(module, "PartElement")
       .def_readonly("nb_owns", &PartElement::nb_owns)
       .def_readonly("nb", &PartElement::nb)
       .def_property_readonly("nb_ghosts", &PartElement::nb_ghosts);

   py::class_<PartInfo>(module, "PartInfo")
       .def_readonly("ncpus", &PartInfo::ncpus)
       .def_readonly("nodes", &PartInfo::nodes)
       .def_readonly("fractures", &PartInfo::fractures)
       .def_readonly("injectors", &PartInfo::injectors)
       .def_readonly("producers", &PartInfo::producers)
       .def_readonly("mswells", &PartInfo::mswells)
       .def_readonly("mswell_nodes", &PartInfo::mswell_nodes);

   module.def("part_info", []() {
      auto res = std::make_unique<PartInfo>();
      MeshSchema_part_info(*res);
      return res;
   });

   //     module.def("SyncPetsc_rowcolnum", [](
   //           py::array_t<int, py::array::c_style>& RowNum,
   //           py::array_t<int, py::array::c_style>& ColNum
   //     ) {
   //           assert(RowNum.ndim()==1);
   //           assert(ColNum.ndim()==1);
   //           SyncPetsc_rowcolnum(
   //                 RowNum.mutable_data(), RowNum.size(),
   //                 ColNum.mutable_data(), ColNum.size()
   //           );
   //     });

   module.def("SyncPetsc_colnum",
              [](py::array_t<int, py::array::c_style>& ColNum) {
                 assert(ColNum.ndim() == 1);
                 SyncPetsc_colnum(ColNum.mutable_data(), ColNum.size());
              });
   module.def("SyncPetscMSWells_colnum",
              [](py::array_t<int, py::array::c_style>& ColNum) {
                 assert(ColNum.ndim() == 1);
                 SyncPetscMSWells_colnum(ColNum.mutable_data(), ColNum.size());
              });

   module.def("SyncPetscMSWells_GetSolMSWell",
              [](py::object V, NewtonIncrements& increments) {
                 auto vec = cast_to_PETSc<Vec>(V);
                 syncpetscmswells_getsolmswell_(&vec, increments.pointers());
              });
}
