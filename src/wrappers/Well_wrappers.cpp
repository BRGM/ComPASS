#include <cassert>
#include <sstream>

#include "Well.h"

#include "Well_wrappers.h"

using namespace ComPASS::Well;

std::ostream& operator<<(std::ostream& os, const Well& well)
{
	auto& geometry = well.geometry;
	os << "Well with radius " << geometry.radius << std::endl;
	os << "and segments:";
	for (auto&& S : geometry.segments) {
		os << std::endl << S[0] << "->" << S[1];
	}
	return os;
}


void add_well_wrappers(py::module& module)
{

	py::class_<Well_geometry>(module, "WellGeometry")
		.def(py::init<>())
		.def_readwrite("radius", &Well_geometry::radius)
		.def("add_segments", [](Well_geometry& instance, py::iterable& segments) {
		auto to_segment = [](py::sequence& s) { 
			return Well_geometry::Segment_type{ { s[0].cast<Well_geometry::Node_id_type>(), s[1].cast<Well_geometry::Node_id_type>() } };
		};
		for (auto&& S : segments) {
			auto seq = S.cast<py::sequence>();
			instance.segments.push_back(to_segment(seq));
		}
	})
		.def("clear_segments", [](Well_geometry& instance) { instance.segments.clear(); });

	py::class_<Well>(module, "Well")
		.def(py::init<>())
		.def_readwrite("geometry", &Well::geometry)
		.def_property_readonly("stopped", &Well::is_stopped)
		.def_property_readonly("injecting", &Well::is_injecting)
		.def_property_readonly("producing", &Well::is_producing)
		.def("__repr__", [](Well& instance) {
		std::ostringstream os;
		os << instance;
		return os.str();
	});

}
