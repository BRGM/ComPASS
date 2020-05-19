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

#include "ArrayWrapper.h"
#include "BlockMatrixWrapper.h"
#include "SyncPetsc_wrappers.h"

namespace py = pybind11;

class LinearSystemBuilder {
   CsrBlockMatrixWrapper JacA;
   DoubleArray2 JacRHS;
   PartInfo myrank_part_info;
   std::vector<PartInfo> part_info;
   std::vector<size_t> rowl_to_rowg;
   std::vector<size_t> coll_to_colg;
   std::size_t n_rowl = 0, n_coll = 0, n_rowg = 0, n_colg = 0;
   std::vector<int> d_nnz;
   std::vector<int> o_nnz;

  public:
   LinearSystemBuilder();
   void compute_nonzeros();
   void compute_ltog();
   auto get_non_zeros() {
      const std::size_t n = static_cast<std::size_t>(n_rowl);
      typedef py::array_t<int, py::array::c_style> nonzeros;
      return py::make_tuple(py::make_tuple(n_rowl, n_rowg),
                            nonzeros{n, d_nnz.data()},
                            nonzeros{n, o_nnz.data()});
   };
   void set_AMPI(py::object pyA);
   void set_RHS(py::object pyRHS);
   std::size_t get_n_rowl() { return n_rowl; };
   std::size_t get_n_coll() { return n_coll; };
   std::size_t get_n_rowg() { return n_rowg; };
   std::size_t get_n_colg() { return n_colg; };
   std::size_t get_d_nnz(size_t i) { return d_nnz[i]; };
   std::size_t get_o_nnz(size_t i) { return o_nnz[i]; };
};
void add_LinearSystem_wrapper(py::module& module);
