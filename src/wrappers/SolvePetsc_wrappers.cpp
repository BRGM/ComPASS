//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <petscvec.h>

#include "ArrayWrapper.h"
#include "PartitioningInformationWrapper.h"
#include "PyBuffer_wrappers.h"
#include "StringWrapper.h"

// Passing PETSc objects between C and Fortran does not rely on iso C binding
// but on opaque pointers and compiler dependant name mangling hence some
// restrictions on the case of the subroutine names. It is supposed that the
// compiler expose the Fortran routine with lower case and a trailing underscore
// (what gcc does). In the above PETSc example a more elaborated strategy is
// used (based on multiple preprocessor defines... !). As this wrapping is bound
// to disappear (one day) we keep it "simple" but not very portable.

extern "C" {
void retrieve_partitioning(PartitioningInformationWrapper&);
void SolvePetsc_SetUp();
int SolvePetsc_KspSolveIterationNumber();
void SolvePetsc_KspSolveIterations(double*, int);
void SolvePetsc_dump_system(const StringWrapper&);
void SolvePetsc_Ksp_configuration(double, int, int);
//     int SolvePetsc_KspSolve();
int compass_petsc_kspsolve_(Vec*);
//     void SolvePetsc_check_solution();
void compass_check_solution_(Vec*);
void retrieve_RowLToRowGBlock(ArrayWrapper&);
void retrieve_ColLToColGBlock(ArrayWrapper&);
void retrieve_RowLToRowG(ArrayWrapper&);
void retrieve_ColLToColG(ArrayWrapper&);
}

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Petsc_caster.h"

void add_SolvePetsc_wrappers(py::module& module) {
   module.def("SolvePetsc_SetUp", &SolvePetsc_SetUp);
   module.def("SolvePetsc_KspSolveIterationNumber",
              &SolvePetsc_KspSolveIterationNumber);
   module.def("SolvePetsc_dump_system", [](py::str basename) {
      SolvePetsc_dump_system(StringWrapper{basename.cast<std::string>()});
   });
   module.def("SolvePetsc_Ksp_configuration", &SolvePetsc_Ksp_configuration);
   module.def("SolvePetsc_Ksp_iterations", []() {
      const auto n = SolvePetsc_KspSolveIterationNumber();
      assert(n >= 0);
      auto res =
          py::array_t<double, py::array::c_style>{static_cast<std::size_t>(n)};
      SolvePetsc_KspSolveIterations(res.mutable_data(), n);
      return res;
   });

   module.def("SolvePetsc_ksp_solve", [](py::object V) {
      auto vec = cast_to_PETSc<Vec>(V);
      return compass_petsc_kspsolve_(&vec);
   });

   module.def("SolvePetsc_check_solution", [](py::object V) {
      auto vec = cast_to_PETSc<Vec>(V);
      compass_check_solution_(&vec);
   });

   module.def(
       "SolvePetsc_RowLToRowGBlock",
       []() { return retrieve_buffer<IntBuffer>(retrieve_RowLToRowGBlock); },
       "Gets the mapping (index set) local row of a physics bloc -> global "
       "row.");

   module.def(
       "SolvePetsc_ColLToColGBlock",
       []() { return retrieve_buffer<IntBuffer>(retrieve_ColLToColGBlock); },
       "Gets the mapping (index set) local column of a physics block -> global "
       "column.");

   module.def(
       "SolvePetsc_RowLToRowG",
       []() { return retrieve_buffer<IntBuffer>(retrieve_RowLToRowG); },
       "Gets the mapping (index set) local row -> global row.");

   module.def(
       "SolvePetsc_ColLToColG",
       []() { return retrieve_buffer<IntBuffer>(retrieve_ColLToColG); },
       "Gets the mapping (index set) local column -> global column.");

   py::class_<PartitioningInformationWrapper>(module, "PartitioningInformation")
       .def("rowl_to_rowg",
            [](const PartitioningInformationWrapper& self) {
               const auto nrl = self.nb_rowl();
               return py::array_t<PartitioningInformationWrapper::integral_type,
                                  py::array::c_style>{nrl, self.rowl_to_rowg};
            })
       .def("coll_to_colg",
            [](const PartitioningInformationWrapper& self) {
               const auto ncl = self.nb_coll();
               return py::array_t<PartitioningInformationWrapper::integral_type,
                                  py::array::c_style>{ncl, self.coll_to_colg};
            })

       .def_readonly("nb_node_own",
                     &PartitioningInformationWrapper::nb_node_own)
       .def_readonly("nb_frac_own",
                     &PartitioningInformationWrapper::nb_frac_own)
       .def_readonly("nb_well_inj_local",
                     &PartitioningInformationWrapper::nb_well_inj_local)
       .def_readonly("nb_well_prod_local",
                     &PartitioningInformationWrapper::nb_well_prod_local)
       .def_readonly("nb_node_local",
                     &PartitioningInformationWrapper::nb_node_local)
       .def_readonly("nb_frac_local",
                     &PartitioningInformationWrapper::nb_frac_local)
       .def_readonly("nb_comp_thermique",
                     &PartitioningInformationWrapper::nb_comp_thermique)
       .def_readonly("nb_well_inj_own",
                     &PartitioningInformationWrapper::nb_well_inj_own)
       .def_readonly("nb_well_prod_own",
                     &PartitioningInformationWrapper::nb_well_prod_own);

   module.def("retrieve_partitioning", []() {
      auto part = std::make_unique<PartitioningInformationWrapper>();
      retrieve_partitioning(*part);
      return part;
   });
}
