//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <petscmat.h>
#include <petscvec.h>

#include "ArrayWrapper.h"
#include "BlockMatrixWrapper.h"
#include "PartitioningInformationWrapper.h"
#include "Petsc_caster.h"
#include "PyBuffer_wrappers.h"
#include "PyXArrayWrapper.h"
#include "StringWrapper.h"
#include "XArrayWrapper.h"
// #include "Jacobian_wrappers.h"

// #include "MeshUtilities.h"

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
void retrieve_jacobian(CsrBlockMatrixWrapper&);
void retrieve_right_hand_side(DoubleArray2&);
void retrieve_nb_nodes_own(XArrayWrapper<int>&);
void retrieve_nb_fractures_own(XArrayWrapper<int>&);
void retrieve_nb_wellinj_own(XArrayWrapper<int>&);
void retrieve_nb_wellprod_own(XArrayWrapper<int>&);
}

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

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

   module.def("set_AMPI_cpp", [](py::object pyA) {
      auto A = cast_to_PETSc<Mat>(pyA);
      // FIXME: we should make a structure that aggregates JacA and part
      CsrBlockMatrixWrapper JacA;
      retrieve_jacobian(JacA);
      PartitioningInformationWrapper part;
      retrieve_partitioning(part);

      const auto n = static_cast<std::size_t>(part.nb_comp_thermique);
      const auto nb_node_own = part.nb_node_own;
      const auto nb_frac_own = part.nb_frac_own;
      const auto nb_node_local = part.nb_node_local;
      const auto nb_frac_local = part.nb_frac_local;
      const auto nb_well_inj_own = part.nb_well_inj_own;
      const auto nb_well_prod_own = part.nb_well_prod_own;

      assert(static_cast<std::size_t>(JacA.block_size) == n);

      std::vector<PetscInt> idxm(
          n, -1);  // number of rows and their global indices
      std::vector<PetscInt> idxn(
          n, -1);  // number of cols and their global indices
      auto block_indices = [n](PetscInt k, std::vector<PetscInt>& idx) {
         for (std::size_t i = 0; i < n; ++i) {
            idx[i] = k;
            ++k;
         }
      };
      auto row_offset = [&JacA](const std::size_t i) {
         return static_cast<std::size_t>(*(JacA.row_offset + i));
      };

      const auto rlg = part.rowl_to_rowg;
      const auto clg = part.coll_to_colg;

      // Rows associated with node own and frac own
      for (std::size_t i = 0; i < nb_node_own + nb_frac_own; ++i) {
         for (auto offset = row_offset(i); offset < row_offset(i + 1);
              ++offset) {
            auto row = *(rlg + i) - 1;  // 0-based in petsc
            auto col =
                *(clg + static_cast<std::size_t>(*(JacA.column + offset) - 1)) -
                1;  // 0-based in petsc
            block_indices(row, idxm);
            block_indices(col, idxn);

            // Col is node or frac, insert JacA.block_data[offset:offset+n*n]
            if (*(JacA.column + offset) - 1 < nb_node_local + nb_frac_local) {
               MatSetValues(A, n, idxm.data(), n, idxn.data(),
                            JacA.block_data + offset * n * n, INSERT_VALUES);
            } else {          // col is wellinj or wellprod, insert
                              // JacA.block_data[offset:offset+n]
               idxn[1] = -1;  // Negative indices are ignored by Petsc.
               MatSetValues(A, n, idxm.data(), n, idxn.data(),
                            JacA.block_data + offset * n * n, INSERT_VALUES);
            };
         }
      }

      // Rows associated with  wellinj own and wellprod own
      // This loop works assuming that the well nodes are put at the end of
      // JacA.block_data
      for (std::size_t i = nb_node_own + nb_frac_own;
           i < nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own;
           ++i) {
         for (auto offset = row_offset(i); offset < row_offset(i + 1);
              ++offset) {
            auto row = *(rlg + i) - 1;
            auto col =
                *(clg + static_cast<std::size_t>(*(JacA.column + offset) - 1)) -
                1;

            block_indices(row, idxm);
            block_indices(col, idxn);

            if (*(JacA.column + offset) - 1 < nb_node_local + nb_frac_local) {
               MatSetValues(A, 1, idxm.data(), n, idxn.data(),
                            JacA.block_data + offset * n * n, INSERT_VALUES);
            } else {
               MatSetValues(A, 1, idxm.data(), 1, idxn.data(),
                            JacA.block_data + offset * n * n, INSERT_VALUES);
            };
         }
      }

      MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
   });
   module.def("get_AMPI_nnz_cpp", []() {
      CsrBlockMatrixWrapper JacA;
      retrieve_jacobian(JacA);
      PartitioningInformationWrapper part;
      retrieve_partitioning(part);

      const auto nb_node_own = part.nb_node_own;
      const auto nb_frac_own = part.nb_frac_own;
      const auto nb_node_local = part.nb_node_local;
      const auto nb_frac_local = part.nb_frac_local;
      const auto nb_comp_thermique = part.nb_comp_thermique;
      const auto nb_well_inj_own = part.nb_well_inj_own;
      const auto nb_well_prod_own = part.nb_well_prod_own;
      const auto nb_well_inj_local = part.nb_well_inj_local;
      const auto nb_well_prod_local = part.nb_well_prod_local;
      XArrayWrapper<int> nb_node_own_ncpus;
      XArrayWrapper<int> nb_frac_own_ncpus;
      XArrayWrapper<int> nb_well_inj_cpus;
      XArrayWrapper<int> nb_well_prod_cpus;
      retrieve_nb_nodes_own(nb_node_own_ncpus);
      retrieve_nb_fractures_own(nb_frac_own_ncpus);
      retrieve_nb_wellinj_own(nb_well_inj_cpus);
      retrieve_nb_wellprod_own(nb_well_prod_cpus);

      auto columns = JacA.column;
      auto row_offset = [&JacA](const std::size_t i) {
         return static_cast<std::size_t>(*(JacA.row_offset + i));
      };

      // local row/col size:  node own and frac own
      auto n_rowl = (nb_node_own + nb_frac_own) * nb_comp_thermique +
                    nb_well_inj_own + nb_well_prod_own;
      auto n_coll = n_rowl;
      // global row/col size: sum of all procs
      size_t n_rowg = 0;
      for (size_t i = 0; i < nb_node_own_ncpus.length; i++) {
         n_rowg +=
             (nb_node_own_ncpus[i] + nb_frac_own_ncpus[i]) * nb_comp_thermique +
             nb_well_inj_cpus[i] + nb_well_prod_cpus[i];
      }
      size_t n_colg = n_rowg;

      // Il y a plein de façon de faire ici...
      // On pourrait en particulier écrire une routine en C++ pur
      // en utilisant srd::vector puis ecrire une fonction d'emballage
      // Ce n'est pas l'option choisi ici...
      // Il y a aussi plein de fa�on d'ecrire l'acces aux donnees
      // auto d_nnz = d_nnz_a.request().ptr
      // peut permttre de garder l'ecriture avec l'arithmetique de pointeur
      // cf. https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
      const std::size_t n = static_cast<std::size_t>(n_rowl);
      auto d_nnz_a = py::array_t<int, py::array::c_style>{n};  // _a for array
      auto d_nnz = d_nnz_a.mutable_unchecked<1>();
      auto o_nnz_a = py::array_t<int, py::array::c_style>{n};
      auto o_nnz = o_nnz_a.mutable_unchecked<1>();
      for (size_t i = 0; i < n; ++i) {
         d_nnz(i) = 0;
         o_nnz(i) = 0;
      }

      // N: Node, F: Frac, WI: well inj, WP: well prod
      // l: local, o: own
      auto nl_fo = nb_node_local + nb_frac_own;
      auto nl_fl = nb_node_local + nb_frac_local;
      auto nl_fl_wio = nb_node_local + nb_frac_local + nb_well_inj_own;
      auto nl_fl_wil = nb_node_local + nb_frac_local + nb_well_inj_local;
      auto nl_fl_wil_wpo =
          nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_own;
      auto nl_fl_wil_wpl = nb_node_local + nb_frac_local + nb_well_inj_local +
                           nb_well_prod_local;

      // rows associated with node own and frac own
      for (size_t i = 0; i < (nb_node_own + nb_frac_own); i++) {
         for (size_t j = row_offset(i); j < row_offset(i + 1); j++) {
            auto lj = columns[j] - 1;

            if (lj < nb_node_own) {  // node own
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  d_nnz(k) += nb_comp_thermique;
               }
            } else if (lj >= nb_node_own && lj < nb_node_local) {  // ghost node
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  o_nnz(k) += nb_comp_thermique;
               }
            } else if (lj >= nb_node_local && lj < nl_fo) {  // frac own
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  d_nnz(k) += nb_comp_thermique;
               }
            } else if (lj >= nl_fo && lj < nl_fl) {  // frac ghost
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  o_nnz(k) += nb_comp_thermique;
               }
            } else if (lj >= nl_fl && lj < nl_fl_wio) {  // well inj own
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  d_nnz(k) += 1;
               }
            } else if (lj >= nl_fl_wio && lj < nl_fl_wil) {  // well inj ghost
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  o_nnz(k) += 1;
               }
            } else if (lj >= nl_fl_wil &&
                       lj < nl_fl_wil_wpo) {  // well prod own
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  d_nnz(k) += 1;
               }
            } else {  // well prod ghost
               for (size_t k = i * nb_comp_thermique;
                    k < (i + 1) * nb_comp_thermique; k++) {
                  o_nnz(k) += 1;
               }
            }
         }
      }

      // rows associated with well inj own and well prod own
      auto start = (nb_node_own + nb_frac_own) * nb_comp_thermique;
      for (size_t i = 0; i < nb_well_inj_own + nb_well_prod_own; i++) {
         auto li = i + nb_node_own + nb_frac_own;

         for (size_t j = row_offset(li); j < row_offset(li + 1); j++) {
            auto lj = columns[j] - 1;

            if (lj < nb_node_own) {  // node own
               d_nnz(start + i) += nb_comp_thermique;
            } else if (nb_node_own <= lj && lj < nb_node_local) {  // node ghost
               o_nnz(start + i) += nb_comp_thermique;

            } else if (nb_node_local <= lj && lj < nl_fo) {  // frac own
               d_nnz(start + i) += nb_comp_thermique;

            } else if (nl_fo <= lj && lj < nl_fl) {  //  frac ghost
               o_nnz(start + i) += nb_comp_thermique;

            } else if (nl_fl <= lj && lj < nl_fl_wio) {  //  well inj own
               d_nnz(start + i) += 1;

            } else if (nl_fl_wio <= lj && lj < nl_fl_wil) {  // well inj ghost
               o_nnz(start + i) += 1;

            } else if (nl_fl_wil <= lj &&
                       lj < nl_fl_wil_wpo) {  //  well prod own
               d_nnz(start + i) += 1;

            } else if (nl_fl_wil_wpo <= lj &&
                       lj < nl_fl_wil_wpl) {  //  well prod ghost
               o_nnz(start + i) += 1;

            } else {
               std::cout << "Error in create Ampi" << std::endl;
            }
         }
      }

      return py::make_tuple(py::make_tuple(n_rowl, n_rowg), d_nnz_a, o_nnz_a);
   });
   module.def("set_RHS_cpp", [](py::object pyRHS) {
      auto RHS = cast_to_PETSc<Vec>(pyRHS);
      DoubleArray2 JacRHS;
      retrieve_right_hand_side(JacRHS);
      PartitioningInformationWrapper part;
      retrieve_partitioning(part);

      const auto nb_node_own = part.nb_node_own;
      const auto nb_frac_own = part.nb_frac_own;
      const auto n = static_cast<std::size_t>(part.nb_comp_thermique);
      const auto nb_well_inj_own = part.nb_well_inj_own;
      const auto nb_well_prod_own = part.nb_well_prod_own;
      const auto rlg = part.rowl_to_rowg;
      int row;
      std::vector<PetscInt> idxm(n, -1);

      auto block_indices = [n](PetscInt k, std::vector<PetscInt>& idx) {
         for (std::size_t i = 0; i < n; ++i) {
            idx[i] = k;
            ++k;
         }
      };

      // Rows associated with node own and frac own
      for (size_t i = 0; i < nb_node_own + nb_frac_own; i++) {
         row = *(rlg + i) - 1;
         block_indices(row, idxm);
         for (size_t j = 0; j < n; j++) {
            VecSetValue(RHS, idxm.data()[j], *(JacRHS.p + n * i + j),
                        INSERT_VALUES);
         }
      }

      // Rows associated with  wellinj own and wellprod own
      // This loop works assuming that the well nodes are put at the end of
      // JacRHS.p
      auto start = nb_node_own + nb_frac_own;
      for (size_t i = 0; i < nb_well_inj_own + nb_well_prod_own; i++) {
         row = *(rlg + i + start) - 1;
         block_indices(row, idxm);
         // VecSetValue takes in a global index, but JacRHS.p is a local
         // variable with block size n
         VecSetValue(RHS, row, *(JacRHS.p + n * (start + i)), INSERT_VALUES);
      }

      VecAssemblyBegin(RHS);
      VecAssemblyEnd(RHS);
   });
}
