#include <cassert>
#include <sstream>

#include "Well.h"

#include "Well_wrappers.h"

using namespace ComPASS::Well;

std::ostream& operator<<(std::ostream& os, const Well& well)
{
	auto& geometry = well.geometry;
	os << "WELL" << std::endl;
	os << "    radius: " << geometry.radius << std::endl;
	os << "    segments: ";
	for (auto&& S : geometry.segments) {
		os << S[0] << "->" << S[1] << ", ";
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
		.def("stop", [](Well& instance) {
		instance.control.status = Stopped_well_status{};
	})
		.def("produce", [](Well& instance) {
		instance.control.status = Production_well_status{};
	})
		.def("inject", [](Well& instance, double T) {
		instance.control.status = Injection_well_status{ T };
	})
		.def_property("operate_on_pressure",
			[](Well& instance) -> py::object {
		py::object result = py::none{};
		if (instance.operates_on_pressure()) {
			const Pressure_operating_conditions * p = boost::get<Pressure_operating_conditions>(&instance.control.operating_conditions);
			assert(p);
			return py::make_tuple(p->pressure, p->flowrate_limit);
		}
		return result;
	},
			[](Well& instance, py::tuple args) {
		instance.control.operating_conditions =
			Pressure_operating_conditions{ args[0].cast<double>(), args[1].cast<double>() };
	})
		.def_property("operate_on_flowrate",
			[](Well& instance) -> py::object {
		py::object result = py::none{};
		if (instance.operates_on_flowrate()) {
			const Flowrate_operating_conditions * p = boost::get<Flowrate_operating_conditions>(&instance.control.operating_conditions);
			assert(p);
			return py::make_tuple(p->flowrate, p->pressure_limit);
		}
		return result;
	},
			[](Well& instance, py::tuple args) {
		instance.control.operating_conditions =
			Flowrate_operating_conditions{ args[0].cast<double>(), args[1].cast<double>() };
	})
		.def("__repr__", [](Well& instance) {
		std::ostringstream os;
		os << instance;
		return os.str();
	});

}
