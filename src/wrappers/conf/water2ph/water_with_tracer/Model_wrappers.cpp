#include "Model_wrappers.h"

void init_model() {}

void finalize_model() {}

void add_model_wrappers(py::module& module) {
   module.def("init_model", &init_model);
   module.def("finalize_model", &finalize_model);
}
