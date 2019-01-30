//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "StringWrapper.h"

// Fortran functions
extern "C"
{
	//void NN_init(const StringWrapper&, const StringWrapper&, const StringWrapper&);
	//void NN_init_warmup_and_read_mesh(const StringWrapper&, const StringWrapper&);
	void NN_init_warmup(const StringWrapper&);
	//void NN_init_read_mesh(const StringWrapper&);
	//void NN_init_build_grid(double, double, double, double, double, double, int, int, int);
	void NN_init_phase2(const StringWrapper&, bool, bool);
	//void NN_main(int, const StringWrapper&);
	// void NN_main_make_timestep(double);
	//void NN_main_output_visu(int, const StringWrapper&);
	void NN_main_summarize_timestep();
	void NN_finalize();
}

#include "NN_wrappers.h"
#include <pybind11/iostream.h>

void add_NN_wrappers(py::module& module)
{

	//module.def("init",
	//	[](const std::string& MeshFile, const std::string& LogFile, const std::string& OutputDir) {
	//	NN_init(MeshFile, LogFile, OutputDir);
	//},
	//	"Initialisation of ComPASS.");

	//module.def("init_warmup_and_read_mesh",
	//[](const std::string& MeshFile, const std::string& LogFile) {
	//	NN_init_warmup_and_read_mesh(MeshFile, LogFile);
	//},
	//"Initialisation of ComPASS up to the point where the global mesh is loaded.  Next logical step is init_phase2.");

	module.def("init_warmup",
	[](const std::string& LogFile) {
		NN_init_warmup(LogFile);
	},
	"Initialisation of ComPASS - warmup phase.");

	//module.def("init_read_mesh",
	//[](const std::string& MeshFile) {
	//	NN_init_read_mesh(MeshFile);
	//},
	//"Initialisation of ComPASS - read the mesh file. Next logical step is init_phase2.");

	//module.def("init_build_grid", &NN_init_build_grid,
	//		   "Initialisation of ComPASS - build a cartesian grid. Next logical step is init_phase2.");

	module.def("init_phase2",
		[](const std::string& OutputDir, bool activate_cpramg, bool activate_direct_solver) {
		NN_init_phase2(OutputDir, activate_cpramg, activate_direct_solver);
	},
		"Initialisation of ComPASS phase 2 : partition and distribute.");

	//module.def("main_loop", [](int TimeIter, const std::string& OutputDir) { NN_main(TimeIter, OutputDir); },
	//	"Main loop of ComPASS.");

	// module.def("make_timestep", &NN_main_make_timestep,
		// py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>());

	// This is transitory to output visualisation files
	//module.def("output_visu", [](int TimeIter, const std::string& OutputDir) { NN_main_output_visu(TimeIter, OutputDir); },
	//	"This function is transitory and is bound to disappear. It is here to output visualisation files through original fortran code.");

	module.def("summarize_timestep", &NN_main_summarize_timestep,
		py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>());

	module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");

}
