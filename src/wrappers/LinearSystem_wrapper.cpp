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

#include "COC.h"
#include "LinearSystemBuilder.h"
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
void retrieve_NumNodebyProc(DofFamilyCOC&);
void retrieve_NumFracbyProc(DofFamilyCOC&);
void retrieve_NumWellProdbyProc(DofFamilyCOC&);
void retrieve_NumWellInjbyProc(DofFamilyCOC&);
void retrieve_NumMSWellNodebyProc(DofFamilyCOC&);
void retrieve_jacobian(CsrBlockMatrixWrapper&);
void retrieve_big_jacobian(CsrBlockMatrixWrapper&);
void retrieve_right_hand_side(DoubleArray2&);
void retrieve_big_right_hand_side(DoubleArray2&);
void MeshSchema_part_info_by_rank(PartInfo&, size_t&);
void MeshSchema_part_info(PartInfo&);
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

/** A function that computes the d_nnz and o_nnz arrays
  d_nnz : array containing the number of nonzeros in the various rows
          of the DIAGONAL portion of the local submatrix (different for each
  row) o_nnz : array containing the number of nonzeros in the various rows of
  the OFF-DIAGONAL portion of the local submatrix (different for each row) See :
  https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateAIJ.html#MatCreateAIJ
*/
void LinearSystemBuilder::compute_nonzeros() {
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
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

   auto columns = JacA.column;
   auto row_offset = [this](const std::size_t i) {
      return static_cast<std::size_t>(*(this->JacA.row_offset + i));
   };

   // local row/col size:  node own and frac own
   n_rowl = (nb_node_own + nb_frac_own) * nb_comp_thermique + nb_well_inj_own +
            nb_well_prod_own + nb_mswell_nodes_own * nb_comp_thermique;
   n_coll = n_rowl;
   // global row/col size: sum of all procs
   n_rowg = 0;
   // rowstarts : proc i possesses rows rowstarts[i] to rowstarts[i+1]-1
   rowstarts.resize(ncpus + 1);
   rowstarts[0] = 0;
   for (size_t i = 0; i < ncpus; i++) {
      n_rowg += (part_info[i].nodes.nb_owns + part_info[i].fractures.nb_owns) *
                    nb_comp_thermique +
                part_info[i].injectors.nb_owns +
                part_info[i].producers.nb_owns +
                part_info[i].mswell_nodes.nb_owns * nb_comp_thermique;

      rowstarts[i + 1] = n_rowg;
   }
   n_colg = n_rowg;

   const std::size_t n = static_cast<std::size_t>(n_rowl);
   d_nnz.clear();
   d_nnz.resize(n, 0);
   o_nnz.clear();
   o_nnz.resize(n, 0);

   // N: Node, F: Frac, WI: well inj, WP: well prod, MWN: mswell-nodes
   // l: local, o: own
   auto nl_fo = nb_node_local + nb_frac_own;
   auto nl_fl = nb_node_local + nb_frac_local;
   auto nl_fl_wio = nb_node_local + nb_frac_local + nb_well_inj_own;
   auto nl_fl_wil = nb_node_local + nb_frac_local + nb_well_inj_local;
   auto nl_fl_wil_wpo =
       nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_own;
   auto nl_fl_wil_wpl =
       nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_local;
   auto nl_fl_wil_wpl_mwno = nb_node_local + nb_frac_local + nb_well_inj_local +
                             nb_well_prod_local + nb_mswell_nodes_own;
   auto nl_fl_wil_wpl_mwnl = nb_node_local + nb_frac_local + nb_well_inj_local +
                             nb_well_prod_local + nb_mswell_nodes_local;

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
         } else if (lj >= nl_fl_wil_wpo &&
                    lj < nl_fl_wil_wpl) {  // well prod ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += 1;
            }
         } else if (lj >= nl_fl_wil_wpl &&
                    lj < nl_fl_wil_wpl_mwno) {  // mswell_node own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[k] += nb_comp_thermique;
            }
         } else if (lj >= nl_fl_wil_wpl_mwno &&
                    lj < nl_fl_wil_wpl_mwnl) {  // mswell_node ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[k] += nb_comp_thermique;
            }
         } else {
            std::cerr << "Error in create Ampi" << std::endl;
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
            std::cerr << "Error in create Ampi" << std::endl;
         }
      }
   }

   // rows associated with mswell-nodes own
   start += nb_well_inj_own + nb_well_prod_own;
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      auto li =
          i + nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own;
      for (size_t j = row_offset(li); j < row_offset(li + 1); j++) {
         auto lj = columns[j] - 1;

         if (lj < nb_node_own) {  //  node own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[start + k] += nb_comp_thermique;
            }
         } else if (nb_node_own <= lj && lj < nb_node_local) {  // node ghost
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[start + k] += nb_comp_thermique;
            }
         } else if (lj >= nl_fl_wil_wpl &&
                    lj < nl_fl_wil_wpl_mwno) {  // mswell node own
            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               d_nnz[start + k] += nb_comp_thermique;
            }
         } else if (lj >= nl_fl_wil_wpl_mwno &&
                    lj < nl_fl_wil_wpl_mwnl) {  // mswell node ghost

            for (size_t k = i * nb_comp_thermique;
                 k < (i + 1) * nb_comp_thermique; k++) {
               o_nnz[start + k] += nb_comp_thermique;
            }
         } else {
            std::cerr << "Error in create Ampi mswellnodes block" << std::endl;
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
void LinearSystemBuilder::compute_ltog() {
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
          part_info[i].injectors.nb_owns + part_info[i].producers.nb_owns +
          (part_info[i].mswell_nodes.nb_owns) * n;
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
   const auto nb_mswell_nodes_own = myrank_part_info.mswell_nodes.nb_owns;
   const auto nb_mswell_nodes_local = myrank_part_info.mswell_nodes.nb;

   rowl_to_rowg.resize(nb_node_own + nb_frac_own + nb_well_inj_own +
                       nb_well_prod_own + nb_mswell_nodes_own);
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

   // ! mswell-nodes
   start += nb_well_inj_own + nb_well_prod_own;
   ndisp += nb_well_inj_own + nb_well_prod_own;
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      rowl_to_rowg[i + start] = ndisp + i * n;
   }

   DofFamilyCOC NBP;
   retrieve_NumNodebyProc(NBP);
   DofFamilyCOC FBP;
   retrieve_NumFracbyProc(FBP);
   DofFamilyCOC WIBP;
   retrieve_NumWellInjbyProc(WIBP);
   DofFamilyCOC WPBP;
   retrieve_NumWellProdbyProc(WPBP);

   DofFamilyCOC MSWN;
   retrieve_NumMSWellNodebyProc(MSWN);

   int local_block_col_offset{0};
   coll_to_colg.resize(nb_node_local + nb_frac_local + nb_well_inj_local +
                       nb_well_prod_local + nb_mswell_nodes_local);

   auto fill_coll_to_colg = [this, &local_block_col_offset, &colstart](
                                DofFamilyCOC dof_family, int dof_size) {
      // A lambda function to automatically fill coll_to_colg, for each unknown
      // type (nodes, fractures, wells, mwells)
      for (size_t i = 0; i < dof_family.size(); i++) {
         auto proc = (dof_family.content_data() + i)->proc;
         auto local_id =
             (dof_family.content_data() + i)
                 ->local_id;  // local_id = dof_family%ids(i)%local_id
         this->coll_to_colg[local_block_col_offset + i] =
             colstart[proc] + (local_id - 1) * dof_size;
      }
   };

   // ! part node local in coll_to_colg
   fill_coll_to_colg(NBP, n);
   local_block_col_offset += nb_node_local;
   for (size_t i = 0; i < ncpus; i++) {
      colstart[i] += n * part_info[i].nodes.nb_owns;
   }

   // // ! part frac in ColLtocolG
   if (nb_frac_local > 0) {
      fill_coll_to_colg(FBP, n);
      local_block_col_offset += nb_frac_local;
      for (size_t i = 0; i < ncpus; i++) {
         colstart[i] += n * part_info[i].fractures.nb_owns;
      }
   }

   // // ! part well inj in ColLtoColG
   if (nb_well_inj_local > 0) {
      fill_coll_to_colg(WIBP, 1);
      local_block_col_offset += nb_well_inj_local;
      for (size_t i = 0; i < ncpus; i++) {
         colstart[i] += part_info[i].injectors.nb_owns;
      }
   }

   // // ! part well prod in ColLtoColG
   if (nb_well_prod_local > 0) {
      fill_coll_to_colg(WPBP, 1);
      local_block_col_offset += nb_well_prod_local;
      for (size_t i = 0; i < ncpus; i++) {
         colstart[i] += part_info[i].producers.nb_owns;
      }
   }

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
void LinearSystemBuilder::set_AMPI(py::object pyA) {
   auto A = cast_to_PETSc<Mat>(pyA);

   const auto n = static_cast<std::size_t>(JacA.block_size);
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_node_local = myrank_part_info.nodes.nb;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_frac_local = myrank_part_info.fractures.nb;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_inj_local = myrank_part_info.injectors.nb;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;
   const auto nb_well_prod_local = myrank_part_info.producers.nb;
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

   // Rows associated with node own and frac own
   for (std::size_t i = 0; i < nb_node_own + nb_frac_own; ++i) {
      for (auto offset = row_offset(i); offset < row_offset(i + 1); ++offset) {
         auto row = rowl_to_rowg[i];
         auto col = coll_to_colg[static_cast<std::size_t>(
             *(JacA.column + offset) - 1)];  // 0-based in petsc
         block_indices(row, idxm);
         block_indices(col, idxn);

         // Col is node,frac or mswell_node, insert
         // JacA.block_data[offset:offset+n*n]
         if ((*(JacA.column + offset) - 1 < nb_node_local + nb_frac_local) ||
             (*(JacA.column + offset) - 1 >= nb_node_local + nb_frac_local +
                                                 nb_well_inj_local +
                                                 nb_well_prod_local)) {
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

   PetscErrorCode ierr;
   // Rows associated with mswell_nodes
   for (std::size_t i =
            nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own;
        i < nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own +
                nb_mswell_nodes_own;
        ++i) {
      for (auto offset = row_offset(i); offset < row_offset(i + 1); ++offset) {
         auto row = rowl_to_rowg[i];
         auto col = coll_to_colg[static_cast<std::size_t>(
             *(JacA.column + offset) - 1)];  // 0-based in petsc
         block_indices(row, idxm);
         block_indices(col, idxn);
         if ((*(JacA.column + offset) - 1 > nb_node_local) &&
             (*(JacA.column + offset) - 1 < nb_node_local + nb_frac_local +
                                                nb_well_inj_local +
                                                nb_well_prod_local)) {
            std::cerr << "Error LinearSystemBuilder::set_AMPI()" << std::endl;

         }

         // Col is node  or mswell_node
         else if (*(JacA.column + offset) - 1 <
                  nb_node_local + nb_frac_local + nb_well_inj_local +
                      nb_well_prod_local + nb_mswell_nodes_local) {
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
void LinearSystemBuilder::set_RHS(py::object pyRHS) {
   auto RHS = cast_to_PETSc<Vec>(pyRHS);

   const auto n = static_cast<std::size_t>(JacRHS.shape[1]);
   const auto nb_node_own = myrank_part_info.nodes.nb_owns;
   const auto nb_frac_own = myrank_part_info.fractures.nb_owns;
   const auto nb_well_inj_own = myrank_part_info.injectors.nb_owns;
   const auto nb_well_prod_own = myrank_part_info.producers.nb_owns;
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

   // Rows associated with mswell_nodes
   start += nb_well_inj_own + nb_well_prod_own;
   for (size_t i = 0; i < nb_mswell_nodes_own; i++) {
      row = rowl_to_rowg[i + start];
      block_indices(row, idxm);
      for (size_t j = 0; j < n; j++) {
         VecSetValue(RHS, idxm.data()[j], *(JacRHS.p + n * (start + i) + j),
                     INSERT_VALUES);
      }
   }

   VecAssemblyBegin(RHS);
   VecAssemblyEnd(RHS);
};

void add_LinearSystem_wrapper(py::module& module) {
   py::class_<LinearSystemBuilder>(module, "LinearSystemBuilder")
       .def(py::init<>())
       .def("get_non_zeros", &LinearSystemBuilder::get_non_zeros)
       .def("get_block_size", &LinearSystemBuilder::get_block_size)
       .def("get_n_wells", &LinearSystemBuilder::get_n_wells)
       .def("get_rowstart", &LinearSystemBuilder::get_rowstart)
       .def("set_AMPI", &LinearSystemBuilder::set_AMPI)
       .def("set_RHS", &LinearSystemBuilder::set_RHS);

   py::class_<CsrBlockMatrixWrapper>(module, "CsrBlockMatrixWrapper")
       .def(
           "offsets",
           [](CsrBlockMatrixWrapper& self) {
              using integral_type = CsrBlockMatrixWrapper::integral_type;
              return py::array_t<integral_type, py::array::c_style>{
                  static_cast<py::ssize_t>(self.nb_rows + 1), self.row_offset};
           },
           py::keep_alive<0, 1>())
       .def(
           "column",
           [](CsrBlockMatrixWrapper& self) {
              using integral_type = CsrBlockMatrixWrapper::integral_type;
              return py::array_t<integral_type, py::array::c_style>{
                  static_cast<py::ssize_t>(self.nb_blocks()), self.column};
           },
           py::keep_alive<0, 1>())
       .def(
           "blocks",
           [](CsrBlockMatrixWrapper& self) {
              using value_type = CsrBlockMatrixWrapper::value_type;
              const auto n = static_cast<py::ssize_t>(self.nb_blocks());
              const auto bs = static_cast<py::ssize_t>(self.block_size);
              constexpr auto u = static_cast<py::ssize_t>(sizeof(value_type));
              return py::array_t<value_type, py::array::c_style>{
                  std::array<py::ssize_t, 3>{n, bs, bs},
                  std::array<py::ssize_t, 3>{bs * bs * u, bs * u, u},
                  self.block_data};
           },
           py::keep_alive<0, 1>())
       .def_readonly("nb_rows", &CsrBlockMatrixWrapper::nb_rows)
       .def_readonly("block_size", &CsrBlockMatrixWrapper::block_size)
       .def("nb_blocks", &CsrBlockMatrixWrapper::nb_blocks);

   module.def("jacobian", []() {
      auto p = std::make_unique<CsrBlockMatrixWrapper>();
      retrieve_jacobian(*p);
      return p;
   });

   module.def("big_jacobian", []() {
      auto p = std::make_unique<CsrBlockMatrixWrapper>();
      retrieve_big_jacobian(*p);
      return p;
   });

   py::class_<DoubleArray2>(module, "DoubleArray2")
       .def(
           "data",  // FIXME: implement buffer interface cf. COC_wrappers.cpp
           [](DoubleArray2& self) {
              const auto n = static_cast<py::ssize_t>(self.shape[0]);
              const auto bs = static_cast<py::ssize_t>(self.shape[1]);
              constexpr auto u = static_cast<py::ssize_t>(sizeof(double));
              return py::array_t<double, py::array::c_style>{
                  std::array<py::ssize_t, 2>{n, bs},
                  std::array<py::ssize_t, 2>{bs * u, u}, self.p};
           },
           py::keep_alive<0, 1>());

   module.def("right_hand_side", []() {
      auto p = std::make_unique<DoubleArray2>();
      retrieve_right_hand_side(*p);
      return p;
   });

   module.def("big_right_hand_side", []() {
      auto p = std::make_unique<DoubleArray2>();
      retrieve_big_right_hand_side(*p);
      return p;
   });
}
