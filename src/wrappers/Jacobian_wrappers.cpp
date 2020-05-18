//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <array>

#include "ArrayWrapper.h"
#include "BlockMatrixWrapper.h"
#include "PyBuffer_wrappers.h"

struct DoubleArray2 {
   double* p;
   std::array<std::size_t, 2> shape;
};

// Fortran functions
extern "C" {
void retrieve_jacobian(CsrBlockMatrixWrapper&);
void retrieve_big_jacobian(CsrBlockMatrixWrapper&);
void retrieve_right_hand_side(DoubleArray2&);
}

#include <pybind11/numpy.h>

#include "Jacobian_wrappers.h"

void add_Jacobian_wrappers(py::module& module) {
   py::class_<CsrBlockMatrixWrapper>(module, "BlockMatrix")
       .def_readonly("nb_rows", &CsrBlockMatrixWrapper::nb_rows)
       .def_readonly("block_size", &CsrBlockMatrixWrapper::block_size)
       .def("blocks",
            [](const CsrBlockMatrixWrapper& self) {
               const auto nb = self.nb_blocks();
               const auto bs = static_cast<std::size_t>(self.block_size);
               return py::array_t<CsrBlockMatrixWrapper::value_type,
                                  py::array::c_style>{{nb, bs, bs},
                                                      self.block_data};
            })
       .def("columns",
            [](const CsrBlockMatrixWrapper& self) {
               const auto nb = self.nb_blocks();
               return py::array_t<CsrBlockMatrixWrapper::integral_type,
                                  py::array::c_style>{nb, self.column};
            })
       // .def("raw_columns", [](const CsrBlockMatrixWrapper& self) {
       //     const auto nb = self.nb_blocks();
       //     return py::array_t<CsrBlockMatrixWrapper::integral_type,
       //     py::array::c_style>{
       //         nb, self.column
       //     };
       // })
       .def("row_offset",
            [](const CsrBlockMatrixWrapper& self) {
               const auto nr = self.nb_rows;
               return py::array_t<CsrBlockMatrixWrapper::integral_type,
                                  py::array::c_style>{
                   static_cast<std::size_t>(nr) + 1, self.row_offset};
            })
       .def("block_indexes", [](const CsrBlockMatrixWrapper& self) {
          const auto nb = self.nb_blocks();
          const auto nr = self.nb_rows;
          typedef CsrBlockMatrixWrapper::integral_type Int;
          auto result = py::array_t<Int, py::array::c_style>{
              {static_cast<std::size_t>(nb), static_cast<std::size_t>(2)}};
          auto pij =
              reinterpret_cast<std::pair<Int, Int>*>(result.mutable_data(0, 0));
          std::size_t k = 0;
          const Int* pj = self.column;
          for (std::size_t i = 0; i != nr; ++i) {
             const auto row_offset =
                 static_cast<std::size_t>(*(self.row_offset + i + 1));
             const auto pj_end = self.column + row_offset;
             for (; pj != pj_end; ++pj) {
                pij->first = static_cast<std::size_t>(i);
                pij->second = static_cast<std::size_t>(*pj);
                ++pij;
             }
          }
          return result;
       });

   module.def("retrieve_jacobian", []() {
      auto A = std::make_unique<CsrBlockMatrixWrapper>();
      retrieve_jacobian(*A);
      return A;
   });

   module.def("retrieve_right_hand_side", []() {
      DoubleArray2 a;
      retrieve_right_hand_side(a);
      return py::array_t<double, py::array::c_style>{a.shape, a.p};
   });
}
