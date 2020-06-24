//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <mpi.h>

// we rely on python interpreter embedding
// cf. http://pybind11.readthedocs.io/en/stable/advanced/embedding.html
#include <pybind11/embed.h>
#include <pybind11/eval.h>

namespace py = pybind11;

int main(int argc, const char* argv[]) {
   MPI_Init(nullptr, nullptr);

   py::scoped_interpreter guard{};
   if (argc != 2) {
      py::print("You must provide exactly one script following the syntax:");
      py::print("\t", argv[0], "path_to_your_script");
      return -1;
   }
   py::object scope = py::module::import("__main__").attr("__dict__");
   scope.attr("setdefault")("__file__", argv[1]);
   py::eval_file(argv[1], scope);

   MPI_Finalize();

   return 0;
}
