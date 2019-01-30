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

#include "NewtonIncrements.h"
#include "CTVector.h"
#include "StringWrapper.h"

// Fortran functions
extern "C"
{
    void IncCV_SaveIncPreviousTimeStep();
    void IncCV_LoadIncPreviousTimeStep();
    void IncCVWells_PressureDrop();
    double IncCVReservoir_NewtonRelax(const NewtonIncrements::Pointers<const double>);
    void IncCV_NewtonIncrement(const NewtonIncrements::Pointers<const double>, const double);
    void DirichletContribution_update();
    void LoisThermoHydro_compute();
    void Flux_DarcyFlux_Cell();
    void Flux_DarcyFlux_Frac();
    void Flux_FourierFlux_Cell();
    void Flux_FourierFlux_Frac();
    void Residu_reset_history();
    void Residu_compute(double);
    double Residu_RelativeNorm_local_closure();
    void Jacobian_ComputeJacSm(double);
    void Jacobian_GetSolCell(NewtonIncrements::Pointers<double>);
    void SolvePetsc_SetUp();
    int SolvePetsc_KspSolveIterationNumber();
    void SolvePetsc_KspSolveIterations(double *, int);
    int SolvePetsc_KspSolve();
    void SolvePetsc_Sync();
    void SolvePetsc_GetSolNodeFracWell(NewtonIncrements::Pointers<double>); 
    void LoisThermoHydro_PrimToSecd(NewtonIncrements::Pointers<double>);
    void NN_flash_all_control_volumes();
    void DefFlashWells_TimeFlash();
    void pass_and_dump_array(double *, std::size_t*);
    void SolvePetsc_dump_system(const StringWrapper&);
    void SolvePetsc_Ksp_configuration(double, int, int);
    void SolvePetsc_check_solution();

}

#include "TimeLoop_wrappers.h"
#include <pybind11/numpy.h>

void add_time_loop_wrappers(py::module& module)
{

    module.def("IncCV_SaveIncPreviousTimeStep", &IncCV_SaveIncPreviousTimeStep);
    module.def("IncCV_LoadIncPreviousTimeStep", &IncCV_LoadIncPreviousTimeStep);
    module.def("IncCVWells_PressureDrop", &IncCVWells_PressureDrop);
    module.def("DirichletContribution_update", &DirichletContribution_update);
    module.def("LoisThermoHydro_compute", &LoisThermoHydro_compute);
    module.def("Flux_DarcyFlux_Cell", &Flux_DarcyFlux_Cell);
    module.def("Flux_DarcyFlux_Frac", &Flux_DarcyFlux_Frac);
    module.def("Flux_FourierFlux_Cell", &Flux_FourierFlux_Cell);
    module.def("Flux_FourierFlux_Frac", &Flux_FourierFlux_Frac);
    module.def("Residu_compute", &Residu_compute);
    module.def("Residu_reset_history", &Residu_reset_history);
    module.def("Residu_RelativeNorm_local_closure", &Residu_RelativeNorm_local_closure);
    module.def("Jacobian_ComputeJacSm", &Jacobian_ComputeJacSm);
    module.def("SolvePetsc_SetUp", &SolvePetsc_SetUp);
    module.def("SolvePetsc_KspSolveIterationNumber", &SolvePetsc_KspSolveIterationNumber);
    module.def("SolvePetsc_KspSolve", &SolvePetsc_KspSolve);
    module.def("SolvePetsc_Sync", &SolvePetsc_Sync);
    module.def("NN_flash_all_control_volumes", &NN_flash_all_control_volumes);
    module.def("DefFlashWells_TimeFlash", &DefFlashWells_TimeFlash);
    module.def("SolvePetsc_check_solution", &SolvePetsc_check_solution);
    module.def("SolvePetsc_dump_system", [](py::str basename) {
        SolvePetsc_dump_system(StringWrapper{ basename.cast<std::string>() });
    });
    module.def("SolvePetsc_Ksp_configuration", &SolvePetsc_Ksp_configuration);
    module.def("SolvePetsc_Ksp_iterations", []() {
        const auto n = SolvePetsc_KspSolveIterationNumber();
        assert(n>=0);
        auto res = py::array_t<double, py::array::c_style>{
            static_cast<std::size_t>(n)
        };
        SolvePetsc_KspSolveIterations(res.mutable_data(), n);
        return res;
    });

    auto as_primary_array = [](std::vector<double>& v) {
        auto res = py::array_t<double, py::array::c_style>{v.size(), v.data()};
        res.attr("shape") = py::make_tuple(-1, NewtonIncrements::npv());
        return res;        
    };
    
    py::class_<NewtonIncrements>(module, "NewtonIncrements")
        .def(py::init())
	.def("init", &NewtonIncrements::init)
	.def("nodes", [as_primary_array](NewtonIncrements& self) {
            return as_primary_array(self.nodes);
		}, py::keep_alive<0, 1>())
	.def("fractures", [as_primary_array](NewtonIncrements& self) {
            return as_primary_array(self.fractures);
		}, py::keep_alive<0, 1>())
	.def("cells", [as_primary_array](NewtonIncrements& self) {
            return as_primary_array(self.cells);
		}, py::keep_alive<0, 1>())
        ;

    py::class_<CTVector>(module, "CTVector")
        .def(py::init())
        .def("as_array", [](CTVector& self) {
			return py::array_t<double, py::array::c_style> {
				self.npv, self.values
			};	
		    },
		py::keep_alive<0, 1>())
    ;

    module.def("IncCVReservoir_NewtonRelax", [](const NewtonIncrements& increments){
	return IncCVReservoir_NewtonRelax(increments.pointers());
      });
    
    module.def("IncCV_NewtonIncrement", [](const NewtonIncrements& increments, const double relaxation){
	IncCV_NewtonIncrement(increments.pointers(), relaxation);
      });
    
    module.def("LoisThermoHydro_PrimToSecd", [](NewtonIncrements& increments){
//     py::print("LT P2S sizes:", nb_primary_variables(), nb_nodes(),
//			    nb_fractures(), nb_cells(), nb_injectors(), nb_producers());
     LoisThermoHydro_PrimToSecd(increments.pointers());
      });
    
    module.def("Jacobian_GetSolCell", [](NewtonIncrements& increments){
 //   py::print("sizes:", nb_primary_variables(), nb_nodes(),
 //       		    nb_fractures(), nb_cells(), nb_injectors(), nb_producers());
     Jacobian_GetSolCell(increments.pointers());
      });
    
    module.def("SolvePetsc_GetSolNodeFracWell", [](NewtonIncrements& increments){
//		    py::print("sizes:", nb_primary_variables(), nb_nodes(),
//				    nb_fractures(), nb_cells(), nb_injectors(), nb_producers());
//		    py::print("      ", increments.nodes.size(), increments.cells.size(), 
//		   increments.injectors.size(), increments.producers.size()); 
		    SolvePetsc_GetSolNodeFracWell(increments.pointers());
       });

}
