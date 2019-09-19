#include <petscmat.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename PetscObject>
auto cast_to_PETSc(py::object obj) {
    auto handle = py::cast<long>(obj.attr("handle"));
    return reinterpret_cast<PetscObject>((void *) handle);
}

PYBIND11_MODULE(mypetsc, module)
{

    module.def("dump", [](py::object mat){
        MatView(cast_to_PETSc<Mat>(mat), PETSC_VIEWER_STDOUT_WORLD);
    });

}
