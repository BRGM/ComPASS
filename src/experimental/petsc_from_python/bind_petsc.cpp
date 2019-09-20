#include <petscmat.h>

// This is based on PETSc example #7:
// https://www.mcs.anl.gov/petsc/petsc-current/src/vec/vec/examples/tutorials/ex7.c
// https://www.mcs.anl.gov/petsc/petsc-current/src/vec/vec/examples/tutorials/ex7f.F
// Passing PETSc objects between C and Fortran does not rely on iso C binding 
// but on opaque pointers and compiler dependant name mangling hence some restrictions
// on the case of the subroutine names.
// It is supposed that the compiler expose the Fortran routine with lower case 
// and a trailing underscore (what gcc does).
// In the above PETSc example a more elaborated strategy is used
// (based on multiple preprocessor defines... !).
// As this wrapping is bound to disappear (one day) we keep it "simple" but not very portable. 

extern "C" {
    // cf. constraints on subroutine names supra
    void dump_from_fortran_(Mat*);
}

#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename PetscObject>
auto cast_to_PETSc(py::object obj) {
    // cython object has a handle attribute that stores the PETSc object adress
    auto handle = py::cast<long>(obj.attr("handle"));
    return reinterpret_cast<PetscObject>((void *) handle);
}

PYBIND11_MODULE(bind_petsc, module)
{

    module.def("dump", [](py::object mat){
        MatView(cast_to_PETSc<Mat>(mat), PETSC_VIEWER_STDOUT_WORLD);
    });
    module.def("dump_from_Fortran", [](py::object obj){
        auto mat = cast_to_PETSc<Mat>(obj);
        dump_from_fortran_(&mat);
    });

}
