//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "ArrayWrapper.h"
#include "StringWrapper.h"

// Fortran functions
extern "C"
{
	//void NN_init(const StringWrapper&, const StringWrapper&, const StringWrapper&);
	//void NN_init_warmup_and_read_mesh(const StringWrapper&, const StringWrapper&);
	void NN_init_warmup(const StringWrapper&);
	//void NN_init_read_mesh(const StringWrapper&);
	//void NN_init_build_grid(double, double, double, double, double, double, int, int, int);
	void NN_init_phase2(bool, bool);
	//void NN_main(int, const StringWrapper&);
	// void NN_main_make_timestep(double);
	//void NN_main_output_visu(int, const StringWrapper&);
	void NN_main_summarize_timestep();
	void NN_finalize();
    void NN_init_phase2_summary();
    void NN_init_phase2_partition(ArrayWrapper&);
}

#include "NN_wrappers.h"

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

void add_NN_wrappers(py::module& module)
{

	module.def("init_warmup",
	[](const std::string& LogFile) {
		NN_init_warmup(LogFile);
	},
	"Initialisation of ComPASS - warmup phase.");

    module.def("init_phase2_summary", &NN_init_phase2_summary,
        "Initial summary of simulation on master proc.");
    module.def("init_phase2_partition", [](py::array_t<int, py::array::c_style> colors) {
        assert(colors.ndim()==1);
        auto wrapper = ArrayWrapper::wrap(colors.mutable_data(), colors.size());
        NN_init_phase2_partition(wrapper);
    },
      "Partition mesh.");
	module.def("init_phase2",
		[](bool activate_cpramg, bool activate_direct_solver) {
		NN_init_phase2(activate_cpramg, activate_direct_solver);
	},
		"Initialisation of ComPASS phase 2 : distribute elements from global to local mesh.");

	module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");

}
