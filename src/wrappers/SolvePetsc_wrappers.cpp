//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "SolvePetsc_wrappers.h"

#include <petscmat.h>
#include <petscvec.h>

#include "Petsc_caster.h"
#include "XArrayWrapper.h"

// Passing PETSc objects between C and Fortran does not rely on iso C binding
// but on opaque pointers and compiler dependant name mangling hence some
// restrictions on the case of the subroutine names. It is assumed that the
// compiler expose the Fortran routine with lower case and a trailing underscore
// (what gcc does). In the above PETSc example a more elaborated strategy is
// used (based on multiple preprocessor defines... !). As this wrapping is bound
// to disappear (one day) we keep it "simple" but not very portable.

extern "C" {
void retrieve_num_node_by_proc(CsrMatrixWrapper&);
void retrieve_num_frac_by_proc(CsrMatrixWrapper&);
void retrieve_num_well_inj_by_proc(CsrMatrixWrapper&);
void retrieve_num_well_prod_by_proc(CsrMatrixWrapper&);
void retrieve_jacobian(CsrBlockMatrixWrapper&);
void retrieve_right_hand_side(DoubleArray2&);
void retrieve_nb_nodes_own(XArrayWrapper<int>&);
void retrieve_nb_fractures_own(XArrayWrapper<int>&);
void retrieve_nb_wellinj_own(XArrayWrapper<int>&);
void retrieve_nb_wellprod_own(XArrayWrapper<int>&);
}

#include <iostream>

LinearSystemBuilder::LinearSystemBuilder() {
   retrieve_jacobian(JacA);
   retrieve_right_hand_side(JacRHS);
   MeshSchema_part_info(myrank_part_info);
   const auto ncpus = myrank_part_info.ncpus;
   part_info.resize(ncpus);
   for (size_t i = 0; i < ncpus; i++) {
      MeshSchema_part_info_by_rank(part_info[i], i);
   }
   compute_nonzeros();
   compute_ltog();
};

void LinearSystemBuilder::compute_nonzeros() {
   // A function that computes the d_nnz and o_nnz arrays
   // d_nnz : array containing the number of nonzeros in the various rows
   //         of the DIAGONAL portion of the local submatrix (different for each
   //         row)
   // o_nnz : array containing the number of nonzeros in the various rows
   //         of the OFF-DIAGONAL portion of the local submatrix (different for
   //         each row)
   //  https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateAIJ.html#MatCreateAIJ
   const auto nb_comp_thermique = JacA.block_size;
   const auto ncpus = myrank_part_info.ncpus;
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_node_local = myrank_part_info.nodes.nb;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_frac_local = myrank_part_info.fractures.nb;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_inj_local = myrank_part_info.injectors.nb;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;
   const auto nb_well_prod_local = myrank_part_info.producers.nb;
   XArrayWrapper<int> nb_node_own_ncpus;
   XArrayWrapper<int> nb_frac_own_ncpus;
   XArrayWrapper<int> nb_well_inj_cpus;
   XArrayWrapper<int> nb_well_prod_cpus;
   retrieve_nb_nodes_own(nb_node_own_ncpus);
   retrieve_nb_fractures_own(nb_frac_own_ncpus);
   retrieve_nb_wellinj_own(nb_well_inj_cpus);
   retrieve_nb_wellprod_own(nb_well_prod_cpus);

   auto columns = JacA.column;
   auto row_offset = [this](const std::size_t i) {
      return static_cast<std::size_t>(*(this->JacA.row_offset + i));
   };

   // local row/col size:  node own and frac own
   n_rowl = (nb_node_own + nb_frac_own) * nb_comp_thermique + nb_well_inj_own +
            nb_well_prod_own;
   n_coll = n_rowl;
   // global row/col size: sum of all procs
   n_rowg = 0;
   for (size_t i = 0; i < nb_node_own_ncpus.length; i++) {
      n_rowg +=
          (nb_node_own_ncpus[i] + nb_frac_own_ncpus[i]) * nb_comp_thermique +
          nb_well_inj_cpus[i] + nb_well_prod_cpus[i];
   }
   n_colg = n_rowg;

   const std::size_t n = static_cast<std::size_t>(n_rowl);
   d_nnz.clear();
   d_nnz.resize(n, 0);
   o_nnz.clear();
   o_nnz.resize(n, 0);

   // N: Node, F: Frac, WI: well inj, WP: well prod
   // l: local, o: own
   auto nl_fo = nb_node_local + nb_frac_own;
   auto nl_fl = nb_node_local + nb_frac_local;
   auto nl_fl_wio = nb_node_local + nb_frac_local + nb_well_inj_own;
   auto nl_fl_wil = nb_node_local + nb_frac_local + nb_well_inj_local;
   auto nl_fl_wil_wpo =
       nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_own;
   auto nl_fl_wil_wpl =
       nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_local;

   // rows associated with node own and frac own
   for (size_t i = 0; i < (nb_node_own + nb_frac_own); i++) {
      for (size_t j = row_offset(i); j < row_offset(i + 1); j++) {
         auto lj = columns[j] - 1;

         if (lj < nb_node_own) {  // node own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nb_node_own && lj < nb_node_local) {  // ghost node
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nb_node_local && lj < nl_fo) {  // frac own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nl_fo && lj < nl_fl) {  // frac ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nl_fl && lj < nl_fl_wio) {  // well inj own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += 1;
            }
         } else if (lj >= nl_fl_wio && lj < nl_fl_wil) {  // well inj ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += 1;
            }
         } else if (lj >= nl_fl_wil && lj < nl_fl_wil_wpo) {  // well prod own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += 1;
            }
         } else {  // well prod ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += 1;
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
            d_nnz[start + i] += nb_comp_thermique;
         } else if (nb_node_own <= lj && lj < nb_node_local) {  // node ghost
            o_nnz[start + i] += nb_comp_thermique;

         } else if (nb_node_local <= lj && lj < nl_fo) {  // frac own
            d_nnz[start + i] += nb_comp_thermique;

         } else if (nl_fo <= lj && lj < nl_fl) {  //  frac ghost
            o_nnz[start + i] += nb_comp_thermique;

         } else if (nl_fl <= lj && lj < nl_fl_wio) {  //  well inj own
            d_nnz[start + i] += 1;

         } else if (nl_fl_wio <= lj && lj < nl_fl_wil) {  // well inj ghost
            o_nnz[start + i] += 1;

         } else if (nl_fl_wil <= lj && lj < nl_fl_wil_wpo) {  //  well prod own
            d_nnz[start + i] += 1;

         } else if (nl_fl_wil_wpo <= lj &&
                    lj < nl_fl_wil_wpl) {  //  well prod ghost
            o_nnz[start + i] += 1;

         } else {
            std::cout << "Error in create Ampi" << std::endl;
         }
      }
   }
};

void LinearSystemBuilder::compute_ltog() {
   // A function that computes the local to global index arrays
   // If i is the local index of a row in the local submatrix
   // then rowl_to_rowg[i] is the global index of that row in the matrix
   // If j is the local index of a column in the local submatrix
   // then coll_to_colg[j] is the global index of that column in the matrix
   const auto ncpus = myrank_part_info.ncpus;
   const auto n = static_cast<std::size_t>(JacA.block_size);
   int rowstart[ncpus];
   int colstart[ncpus];
   rowstart[0] = 0;
   colstart[0] = rowstart[0];
   size_t i;
   for (i = 0; i < ncpus - 1; i++) {
      rowstart[i + 1] =
          rowstart[i] +
          (part_info[i].nodes.nb_owns + part_info[i].fractures.nb_owns) * n +
          part_info[i].injectors.nb_owns + part_info[i].producers.nb_owns;
      colstart[i + 1] = rowstart[i + 1];
   }

   const auto comm_rank = myrank_part_info.rank;
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_node_local = myrank_part_info.nodes.nb;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_frac_local = myrank_part_info.fractures.nb;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_inj_local = myrank_part_info.injectors.nb;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;
   const auto nb_well_prod_local = myrank_part_info.producers.nb;

   auto nblock_rowl = static_cast<std::size_t>(
       nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own);
   rowl_to_rowg.resize(nblock_rowl);
   // ! node/frac
   for (size_t i = 0; i < nb_node_own + nb_frac_own; i++) {
      rowl_to_rowg[i] = rowstart[comm_rank] + i * n;  // C indexing
   }

   // ! well
   auto start = nb_node_own + nb_frac_own;
   auto ndisp = rowstart[comm_rank] + (nb_node_own + nb_frac_own) * n;

   for (size_t i = 0; i < nb_well_inj_own + nb_well_prod_own; i++) {
      rowl_to_rowg[i + start] = ndisp + i;
   }

   coll_to_colg.resize(nb_node_local + nb_frac_local + nb_well_inj_local +
                       nb_well_prod_local);
   CsrMatrixWrapper NBP;
   retrieve_num_node_by_proc(NBP);
   auto nb = NBP.nb_rows;
   CsrMatrixWrapper FBP;
   retrieve_num_frac_by_proc(FBP);
   CsrMatrixWrapper WIBP;
   retrieve_num_well_inj_by_proc(WIBP);
   CsrMatrixWrapper WPBP;
   retrieve_num_well_prod_by_proc(WPBP);

   // ! part node local in ColLtocolG
   // Nb : nb_rows, data : Val, column : Num, row_offset : Pt
   for (size_t i = 0; i < NBP.nb_rows; i++) {
      for (size_t j = NBP.row_offset[i]; j < NBP.row_offset[i + 1]; j++) {
         // ! j is in proc ipc
         auto ipc = NBP.data[j];
         coll_to_colg[j] = (NBP.column[j] - 1) * n + colstart[ipc];
      }
   }

   // ! part frac in ColLtocolG
   for (size_t i = 0; i < FBP.nb_rows; i++) {
      for (size_t j = FBP.row_offset[i]; j < FBP.row_offset[i + 1]; j++) {
         //       ! j is in proc ipc
         auto ipc = FBP.data[j];
         coll_to_colg[j + nb_node_local] = (FBP.column[j] - 1) * n +
                                           colstart[ipc] +
                                           part_info[ipc].nodes.nb_owns * n;
      }
   }

   // ! part well inj in ColLtoColG
   start = nb_node_local + nb_frac_local;
   for (size_t i = 0; i < WIBP.nb_rows; i++) {
      for (size_t j = WIBP.row_offset[i]; j < WIBP.row_offset[i + 1]; j++) {
         //       ! j is in proc ipc
         auto ipc = WIBP.data[j];
         coll_to_colg[j + start] =
             WIBP.column[j] - 1 + colstart[ipc] +
             (part_info[ipc].nodes.nb_owns + part_info[ipc].fractures.nb_owns) *
                 n;
         // std::cout << j + start << "  " << coll_to_colg[j + start] <<
         // std::endl;
      }
   }

   // ! part well prod in ColLtoColG
   start = nb_node_local + nb_frac_local + nb_well_inj_local;
   for (size_t i = 0; i < WPBP.nb_rows; i++) {
      for (size_t j = WPBP.row_offset[i]; j < WPBP.row_offset[i + 1]; j++) {
         //  j is in proc ipc
         auto ipc = WPBP.data[j];
         coll_to_colg[j + start] =
             WPBP.column[j] - 1 + colstart[ipc] +
             (part_info[ipc].nodes.nb_owns + part_info[ipc].fractures.nb_owns) *
                 n +
             part_info[ipc].injectors.nb_owns;
         // std::cout << j + start << "  " << coll_to_colg[j + start] <<
         // std::endl;
      }
   }
};

void LinearSystemBuilder::set_AMPI(py::object pyA) {
   // A function that sets the nonzeros entries in the Petsc matrix
   // using the JacA.block_data array computed in the Fortran core.
   // This is done in the C++ layer for increased performance
   auto A = cast_to_PETSc<Mat>(pyA);

   const auto n = static_cast<std::size_t>(JacA.block_size);
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_node_local = myrank_part_info.nodes.nb;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_frac_local = myrank_part_info.fractures.nb;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;

   std::vector<PetscInt> idxm(n,
                              -1);  // number of rows and their global indices
   std::vector<PetscInt> idxn(n,
                              -1);  // number of cols and their global indices
   auto block_indices = [n](PetscInt k, std::vector<PetscInt>& idx) {
      for (std::size_t i = 0; i < n; ++i) {
         idx[i] = k;
         ++k;
      }
   };
   auto row_offset = [this](const std::size_t i) {
      return static_cast<std::size_t>(*(this->JacA.row_offset + i));
   };

   // Rows associated with node own and frac own
   for (std::size_t i = 0; i < nb_node_own + nb_frac_own; ++i) {
      for (auto offset = row_offset(i); offset < row_offset(i + 1); ++offset) {
         auto row = rowl_to_rowg[i];
         auto col = coll_to_colg[static_cast<std::size_t>(
             *(JacA.column + offset) - 1)];  // 0-based in petsc
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
      for (auto offset = row_offset(i); offset < row_offset(i + 1); ++offset) {
         auto row = rowl_to_rowg[i];
         auto col = coll_to_colg[static_cast<std::size_t>(
             *(JacA.column + offset) - 1)];

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
};

void LinearSystemBuilder::set_RHS(py::object pyRHS) {
   // A function that sets the nonzeros entries in the Petsc vector RHS
   // using the JacRHS.p array computed in the Fortran core.
   // This is done in the C++ layer for increased performance
   auto RHS = cast_to_PETSc<Vec>(pyRHS);

   const auto n = static_cast<std::size_t>(JacRHS.shape[1]);
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;

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
      row = rowl_to_rowg[i];
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
      row = rowl_to_rowg[i + start];
      block_indices(row, idxm);
      // VecSetValue takes in a global index, but JacRHS.p is a local variable
      // with block size n
      VecSetValue(RHS, row, *(JacRHS.p + n * (start + i)), INSERT_VALUES);
   }

   VecAssemblyBegin(RHS);
   VecAssemblyEnd(RHS);
};

void add_SolvePetsc_wrappers(py::module& module) {
   py::class_<LinearSystemBuilder>(module, "LinearSystemBuilder")
       .def(py::init<>())
       .def("get_non_zeros", &LinearSystemBuilder::get_non_zeros)
       .def("set_AMPI", &LinearSystemBuilder::set_AMPI)
       .def("set_RHS", &LinearSystemBuilder::set_RHS);
}
