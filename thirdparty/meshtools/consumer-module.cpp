#include <vector>

#include <pybind11/pybind11.h>
//// mandatory to have optional_caster and variant_caster
//#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "meshtools.h"

namespace py = pybind11;
namespace MT = MeshTools;


// The following makes all int vectors opaque
// (cf. PYBIND11_MAKE_OPAQUE)
namespace pybind11 {
	namespace detail {
		template<>
		class type_caster<std::vector<MT::NodeId>> :
			public type_caster_base<std::vector<int>>
		{};
	}
}

PYBIND11_MODULE(consumer, module)
{

	module.doc() = "pybind11 mesh objects consumer test module";

	py::bind_vector<std::vector<MT::NodeId>>(module, "IdVector");

	module.def("consume_vector", [](const std::vector<int>& v) {
		py::print("Vector:");
		for (auto&& x : v) {
			py::print("", x);
		}
	});

}

