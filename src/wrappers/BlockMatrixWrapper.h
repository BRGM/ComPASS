//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <cassert>
#include <cstddef>

struct DoubleArray2 {
   double* p;
   std::array<std::size_t, 2> shape;
};

struct CsrBlockMatrixWrapper {
   typedef double value_type;
   typedef int integral_type;
   const integral_type* row_offset;
   const integral_type* column;
   const value_type* block_data;
   const integral_type nb_rows;
   const integral_type block_size;
   CsrBlockMatrixWrapper()
       : row_offset{nullptr},
         column{nullptr},
         block_data{nullptr},
         nb_rows{0},
         block_size{0} {}
   auto nb_blocks() const {
      assert(row_offset);
      return static_cast<std::size_t>(*(row_offset + nb_rows));
   }
};

struct CsrMatrixWrapper {
   typedef int value_type;
   typedef int integral_type;
   const integral_type* row_offset;
   const integral_type* column;
   const value_type* data;
   const integral_type nb_rows;
   CsrMatrixWrapper()
       : row_offset{nullptr}, column{nullptr}, data{nullptr}, nb_rows{0} {}
   auto nb_blocks() const {
      assert(row_offset);
      return static_cast<std::size_t>(*(row_offset + nb_rows));
   }
};
