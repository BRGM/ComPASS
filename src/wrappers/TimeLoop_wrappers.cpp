//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

//#include "ArrayWrapper.h"
//#include "XArrayWrapper.h"
//#include "PyXArrayWrapper.h"
//#include "PyArrayWrapper.h"
//#include "PyBuffer_wrappers.h"
//#include "MeshUtilities.h"
#include <cstddef>

// Fortran functions
extern "C"
{
    void IncCV_SaveIncPreviousTimeStep();
    void IncCV_LoadIncPreviousTimeStep();
    void IncCVWells_PressureDrop();
    void IncCVReservoir_NewtonRelax(); // NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax = double& -> function
    void IncCV_NewtonIncrement(); // everything intent in NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, NewtonIncreWellInj, NewtonIncreWellProd, relax
    void DirichletContribution_update();
    void LoisThermoHydro_compute();
    void Flux_DarcyFlux_Cell();
    void Flux_DarcyFlux_Frac();
    void Flux_FourierFlux_Cell();
    void Flux_FourierFlux_Frac();
    void Residu_compute(double, int);
    void Residu_RelativeNorm(int, double, double&, double&, double&);
    void Jacobian_ComputeJacSm(int);
    void Jacobian_GetSolCell(); // NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell
    void SolvePetsc_SetUp();
    void SolvePetsc_KspSolve(); // NkspIter, kspHistoryOutput
    void SolvePetsc_Sync();
    void SolvePetsc_GetSolNodeFracWell(); // NewtonIncreNode, NewtonIncreFrac, NewtonIncreWellInj, NewtonIncreWellProd
    void LoisThermoHydro_PrimToSecd(); // vnode, vfrac, vcell
    void NN_flash_all_control_volumes();
    void DefFlashWells_TimeFlash();
    void pass_and_dump_array(double *, std::size_t*);
    //void pass_and_dump_dim_array(double *, int, std::size_t*);
}

#include "TimeLoop_wrappers.h"

void add_time_loop_wrappers(py::module& module)
{

    module.def("IncCV_SaveIncPreviousTimeStep", &IncCV_SaveIncPreviousTimeStep);
    module.def("IncCV_LoadIncPreviousTimeStep", &IncCV_LoadIncPreviousTimeStep);
    module.def("IncCVWells_PressureDrop", &IncCVWells_PressureDrop);
    module.def("IncCVReservoir_NewtonRelax", &IncCVReservoir_NewtonRelax);
    module.def("IncCV_NewtonIncrement", &IncCV_NewtonIncrement);
    module.def("DirichletContribution_update", &DirichletContribution_update);
    module.def("LoisThermoHydro_compute", &LoisThermoHydro_compute);
    module.def("Flux_DarcyFlux_Cell", &Flux_DarcyFlux_Cell);
    module.def("Flux_DarcyFlux_Frac", &Flux_DarcyFlux_Frac);
    module.def("Flux_FourierFlux_Cell", &Flux_FourierFlux_Cell);
    module.def("Flux_FourierFlux_Frac", &Flux_FourierFlux_Frac);
    module.def("Residu_compute", &Residu_compute);
    module.def("Residu_RelativeNorm", &Residu_RelativeNorm);
    module.def("Jacobian_ComputeJacSm", &Jacobian_ComputeJacSm);
    module.def("Jacobian_GetSolCell", &Jacobian_GetSolCell);
    module.def("SolvePetsc_SetUp", &SolvePetsc_SetUp);
    module.def("SolvePetsc_KspSolve", &SolvePetsc_KspSolve);
    module.def("SolvePetsc_Sync", &SolvePetsc_Sync);
    module.def("SolvePetsc_GetSolNodeFracWell", &SolvePetsc_GetSolNodeFracWell);
    module.def("LoisThermoHydro_PrimToSecd", &LoisThermoHydro_PrimToSecd);
    module.def("NN_flash_all_control_volumes", &NN_flash_all_control_volumes);
    module.def("DefFlashWells_TimeFlash", &DefFlashWells_TimeFlash);

    module.def("test_pass_and_dump_array", [] {
        std::vector<double> v;
        for (int i = 0; i < 8; ++i) v.push_back(i);
        std::vector<std::size_t> shape;
        shape.push_back(4);
        shape.push_back(2);
        pass_and_dump_array(v.data(), shape.data());
        //pass_and_dump_dim_array(v.data(), 2, shape.data());
    });

}
