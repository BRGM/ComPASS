// Fortran functions
extern "C"
{
	double get_current_time();
	double get_delta_t();
	double get_final_time();
	void set_final_time(double);
}

#include "GlobalVariables_wrappers.h"

void add_global_variables_wrappers(py::module& module)
{

	module.def("get_current_time", &get_current_time);
	module.def("get_timestep", &get_delta_t);
	module.def("get_final_time", &get_final_time);
	module.def("set_final_time", &set_final_time);

}