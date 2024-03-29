//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <limits>

#include "BlockMatrixWrapper.h"
#include "SyncPetsc_wrappers.h"

namespace py = pybind11;

class LinearSystemBuilderMSWells {
   CsrBlockMatrixWrapper JacA;
   DoubleArray2 JacRHS;
   PartInfo myrank_part_info;
   std::vector<PartInfo> part_info;
   std::vector<size_t> rowl_to_rowg;
   std::vector<size_t> coll_to_colg;
   std::size_t n_rowl = 0, n_coll = 0, n_rowg = 0, n_colg = 0;
   std::vector<int> d_nnz;
   std::vector<int> o_nnz;
   std::vector<size_t> rowstarts;

  public:
   LinearSystemBuilderMSWells();
   void compute_nonzeros();
   void compute_ltog();
   auto get_non_zeros() {
      const std::size_t n = static_cast<std::size_t>(n_rowl);
      typedef py::array_t<int, py::array::c_style> nonzeros;
      assert(n < std::numeric_limits<py::ssize_t>::max());
      return py::make_tuple(
          py::make_tuple(n_rowl, n_rowg),
          nonzeros{static_cast<py::ssize_t>(n), d_nnz.data()},
          nonzeros{static_cast<py::ssize_t>(n), o_nnz.data()});
   };
   void set_AMPI(py::object pyA);
   void set_RHS(py::object pyRHS);
   std::size_t get_block_size() { return JacA.block_size; };
   std::size_t get_n_rowl() { return n_rowl; };
   std::size_t get_n_coll() { return n_coll; };
   std::size_t get_n_rowg() { return n_rowg; };
   std::size_t get_n_colg() { return n_colg; };
   std::size_t get_d_nnz(size_t i) { return d_nnz[i]; };
   std::size_t get_o_nnz(size_t i) { return o_nnz[i]; };
   //   std::size_t get_n_wells() {
   //      return myrank_part_info.producers.nb_owns +
   //             myrank_part_info.injectors.nb_owns;
   //   };

   std::size_t get_rowstart(size_t proc) { return rowstarts[proc]; };
};
