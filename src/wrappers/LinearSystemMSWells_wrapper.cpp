//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include "COC.h"
#include "LinearSystemBuilderMSWells.h"
#include "Petsc_caster.h"
#include "XArrayWrapper.h"

// Passing PETSc objects between C and Fortran does not rely on iso C binding
// but on opaque pointers and compiler dependant name mangling hence some
// restrictions on the case of the subroutine names. It is assumed that the
// compiler expose the Fortran routine with lower case and a trailing underscore
// (what gcc does). In the above PETSc example a more elaborated strategy is
// used (based on multiple preprocessor defines... !). As this wrapping is bound
// to disappear (one day) we keep it "simple" but not very portable.

struct DOFId {
   int proc;
   int local_id;
};

typedef GenericCOC<DOFId> DofFamilyCOC;

// Fortran functions

extern "C" {
void retrieve_jacobian_mswells(CsrBlockMatrixWrapper&);
void retrieve_right_hand_side_mswells(DoubleArray2&);
void MeshSchema_part_info_by_rank(PartInfo&, size_t&);
void MeshSchema_part_info(PartInfo&);
void retrieve_NumMSWellNodebyProc(DofFamilyCOC&);
}

#include <iostream>

LinearSystemBuilderMSWells::LinearSystemBuilderMSWells() {
   retrieve_jacobian_mswells(JacA);
   retrieve_right_hand_side_mswells(JacRHS);
   MeshSchema_part_info(myrank_part_info);
   const auto ncpus = myrank_part_info.ncpus;
   part_info.resize(ncpus);
   for (size_t i = 0; i < ncpus; i++) {
      MeshSchema_part_info_by_rank(part_info[i], i);
   }
   compute_nonzeros();
   compute_ltog();
};

/** A function that computes the d_nnz and o_nnz arrays
  d_nnz : array containing the number of nonzeros in the various rows
          of the DIAGONAL portion of the local submatrix (different for each
  row) o_nnz : array containing the number of nonzeros in the various rows of
  the OFF-DIAGONAL portion of the local submatrix (different for each row) See :
  https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateAIJ.html#MatCreateAIJ
*/
void LinearSystemBuilderMSWells::compute_nonzeros() {
   const auto nb_comp_thermique = JacA.block_size;
   const auto ncpus = myrank_part_info.ncpus;
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

   auto columns = JacA.column;
   auto row_offset = [this](const std::size_t i) {
      return static_cast<std::size_t>(*(this->JacA.row_offset + i));
   };

   // local row/col size:  mswell-nodes own
   n_rowl = nb_mswell_nodes_own * nb_comp_thermique;
   n_coll = n_rowl;
   // global row/col size: sum of all procs
   n_rowg = 0;
   for (size_t i = 0; i < ncpus; i++) {
      n_rowg += part_info[i].mswell_nodes.nb_owns * nb_comp_thermique;
   }
   n_colg = n_rowg;

   const std::size_t n = static_cast<std::size_t>(n_rowl);
   d_nnz.clear();
   d_nnz.resize(n, 0);
   o_nnz.clear();
   o_nnz.resize(n, 0);

   // N: Node, F: Frac, WI: well inj, WP: well prod, MWN: mswell-nodes
   // l: local, o: own
   auto nl_mwno = nb_mswell_nodes_own;
   auto nl_mwnl = nb_mswell_nodes_local;

   // rows associated with mswell-nodes own
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      for (size_t j = row_offset(i); j < row_offset(i + 1); j++) {
         auto lj = columns[j] - 1;

         if (lj < nb_mswell_nodes_own) {  // node own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nb_mswell_nodes_own &&
                    lj < nb_mswell_nodes_local) {  // ghost node
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += nb_comp_thermique;
            }
         } else {
            std::cerr << "Error in create Ampi" << std::endl;
         }
      }
   }
};

/** A function that computes the local to global index arrays
    If i is the local index of a row in the local submatrix
    then rowl_to_rowg[i] is the global index of that row in the matrix
    If j is the local index of a column in the local submatrix
    then coll_to_colg[j] is the global index of that column in the matrix
*/
void LinearSystemBuilderMSWells::compute_ltog() {
   const auto ncpus = myrank_part_info.ncpus;
   const auto n = static_cast<std::size_t>(JacA.block_size);
   int rowstart[ncpus];
   int colstart[ncpus];
   rowstart[0] = 0;
   colstart[0] = rowstart[0];
   size_t i;

   for (i = 0; i < ncpus - 1; i++) {
      rowstart[i + 1] = rowstart[i] + (part_info[i].mswell_nodes.nb_owns) * n;
      colstart[i + 1] = rowstart[i + 1];
   }

   const auto comm_rank = myrank_part_info.rank;
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

   rowl_to_rowg.resize(nb_mswell_nodes_own);
   // ! mswell-nodes
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      rowl_to_rowg[i] = rowstart[comm_rank] + i * n;  // C indexing
   }

   DofFamilyCOC MSWN;
   retrieve_NumMSWellNodebyProc(MSWN);

   int local_block_col_offset{0};
   coll_to_colg.resize(nb_mswell_nodes_local);

   auto fill_coll_to_colg = [this, &local_block_col_offset, &colstart](
                                DofFamilyCOC dof_family, int dof_size) {
      // A lambda function to automatically fill coll_to_colg, for each unknown
      // type (nodes, fractures, wells)
      for (size_t i = 0; i < dof_family.size(); i++) {
         auto proc = (dof_family.content_data() + i)->proc;
         auto local_id =
             (dof_family.content_data() + i)
                 ->local_id;  // local_id = dof_family%ids(i)%local_id
         this->coll_to_colg[local_block_col_offset + i] =
             colstart[proc] + (local_id - 1) * dof_size;
      }
   };

   // ! part mswell_nodes local in coll_to_colg
   if (nb_mswell_nodes_local > 0) {
      fill_coll_to_colg(MSWN, n);
      local_block_col_offset += nb_mswell_nodes_local;
      for (size_t i = 0; i < ncpus; i++) {
         colstart[i] += n * part_info[i].mswell_nodes.nb_owns;
      }
   }
};

/** A function that sets the nonzeros entries in the Petsc matrix
   using the JacA.block_data array computed in the Fortran core.
   This is done in the C++ layer for increased performance
*/

void LinearSystemBuilderMSWells::set_AMPI(py::object pyA) {
   auto A = cast_to_PETSc<Mat>(pyA);

   const auto n = static_cast<std::size_t>(JacA.block_size);
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

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
   PetscErrorCode ierr;

   // Rows associated with mswell_nodes
   for (std::size_t i = 0; i < nb_mswell_nodes_own; ++i) {
      for (auto offset = row_offset(i); offset < row_offset(i + 1); ++offset) {
         auto row = rowl_to_rowg[i];
         auto col = coll_to_colg[static_cast<std::size_t>(
             *(JacA.column + offset) - 1)];  // 0-based in petsc
         block_indices(row, idxm);
         block_indices(col, idxn);

         // Col is node or frac, insert JacA.block_data[offset:offset+n*n]
         if (*(JacA.column + offset) - 1 < nb_mswell_nodes_local) {
            ierr =
                MatSetValues(A, n, idxm.data(), n, idxn.data(),
                             JacA.block_data + offset * n * n, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
         }
      }
   }

   ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
   CHKERRABORT(PETSC_COMM_WORLD, ierr);
   ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
   CHKERRABORT(PETSC_COMM_WORLD, ierr);
};

/** A function that sets the nonzeros entries in the Petsc vector RHS
    using the JacRHS.p array computed in the Fortran core.
    This is done in the C++ layer for increased performance
*/
void LinearSystemBuilderMSWells::set_RHS(py::object pyRHS) {
   auto RHS = cast_to_PETSc<Vec>(pyRHS);

   const auto n = static_cast<std::size_t>(JacRHS.shape[1]);
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

   int row;
   std::vector<PetscInt> idxm(n, -1);

   auto block_indices = [n](PetscInt k, std::vector<PetscInt>& idx) {
      for (std::size_t i = 0; i < n; ++i) {
         idx[i] = k;
         ++k;
      }
   };

   // Rows associated with mswell_nodes
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      row = rowl_to_rowg[i];
      block_indices(row, idxm);
      for (size_t j = 0; j < n; j++) {
         VecSetValue(RHS, idxm.data()[j], *(JacRHS.p + n * i + j),
                     INSERT_VALUES);
      }
   }

   VecAssemblyBegin(RHS);
   VecAssemblyEnd(RHS);
};

void add_LinearSystemMSWells_wrapper(py::module& module) {
   py::class_<LinearSystemBuilderMSWells>(module, "LinearSystemBuilderMSWells")
       .def(py::init<>())
       .def("get_non_zeros", &LinearSystemBuilderMSWells::get_non_zeros)
       .def("get_block_size", &LinearSystemBuilderMSWells::get_block_size)
       //.def("get_n_wells", &LinearSystemBuilderMSWells::get_n_wells)
       .def("get_rowstart", &LinearSystemBuilderMSWells::get_rowstart)
       .def("set_AMPI", &LinearSystemBuilderMSWells::set_AMPI)
       .def("set_RHS", &LinearSystemBuilderMSWells::set_RHS);
}
