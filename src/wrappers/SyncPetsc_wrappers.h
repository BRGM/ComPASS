//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

struct MatrixSize {
   typedef std::size_t size_type;
   size_type nb_rows;
   size_type nb_cols;
};

struct PartElement {
   typedef std::size_t size_type;
   size_type nb_owns;
   size_type nb;
   size_type nb_ghosts() const {
      assert(nb_owns <= nb);
      return nb - nb_owns;
   }
};

struct PartInfo {
   typedef std::size_t size_type;
   size_type ncpus;
   size_type rank;
   PartElement nodes;
   PartElement fractures;
   PartElement injectors;
   PartElement producers;
};

extern "C" {
void MeshSchema_part_info_by_rank(PartInfo&, size_t&);
void MeshSchema_part_info(PartInfo&);
}

void add_SyncPetsc_wrappers(py::module& module);
