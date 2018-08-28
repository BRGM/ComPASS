//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <cassert>
#include <sstream>
#include <functional>
#include <iostream>

#include "ArrayWrapper.h"
#include "COC.h"
#include "Well.h"

#include "Well_wrappers.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;

/** Common data structure shared by injectors and producers. */
struct Fortran_well_data
{
    typedef std::array<double, NC> Component_vector;
    char operating_code;
    double radius;
    double maximum_pressure;
    double minimum_pressure;
    double imposed_flowrate;
    Component_vector injection_composition;
    double injection_temperature;
};

// Fortran functions
extern "C"
{
	void Well_allocate_well_geometries(COC&, COC&);
	void Well_set_wells_data(ArrayWrapper&, ArrayWrapper&);
    Fortran_well_data * get_injectors_data();
    int nb_injectors();
    Fortran_well_data * get_producers_data();
    int nb_producers();
}

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

struct WellGeometryInfo {
	std::vector<int> offsets;
	std::vector<int> edges;
	void shrink_to_fit() {
		offsets.shrink_to_fit();
		edges.shrink_to_fit();
	}
};

auto extract_geometry_info(py::list wells, std::function<bool(const Well&)> select)
{
	auto geometry_info = WellGeometryInfo{};
	auto& offsets = geometry_info.offsets;
	auto& edges = geometry_info.edges;
	assert(edges.empty());
	offsets.push_back(0);
	for (auto&& p : wells) {
		auto well = p.cast<const Well&>();
		if (select(well)) {
			for (auto&& edge : well.geometry.segments) {
				edges.push_back(edge[0]);
				edges.push_back(edge[1]);
			}
			offsets.push_back(edges.size());
		}
	}
	geometry_info.shrink_to_fit();
	return geometry_info;
}

void set_well_geometries(py::list wells)
{
	auto producers_geometries = extract_geometry_info(wells, [](const Well& well) { return well.is_producing(); });
	auto injectors_geometries = extract_geometry_info(wells, [](const Well& well) { return well.is_injecting(); });
	auto producers = COC::wrap(producers_geometries.offsets, producers_geometries.edges);
	auto injectors = COC::wrap(injectors_geometries.offsets, injectors_geometries.edges);
	Well_allocate_well_geometries(producers, injectors);
}

struct Producer_data_info
{
	double radius;
	double minimum_pressure;
	double imposed_flowrate;
	char operating_code;  // 'p' for pressure mode; 'f' for flowrate mode
	Producer_data_info(const Well& well) {
		assert(well.is_producing());
		radius = well.geometry.radius;
		if (well.operates_on_pressure()) {	
			operating_code = 'p';
			assert(well.control.operating_conditions.is<Pressure_operating_conditions>());
			auto operating_conditions = well.control.operating_conditions.get<Pressure_operating_conditions>();
			minimum_pressure = operating_conditions.pressure;
			imposed_flowrate = operating_conditions.flowrate_limit;
		}
		else {
			assert(well.operates_on_flowrate());
			operating_code = 'f';
			assert(well.control.operating_conditions.is<Flowrate_operating_conditions>());
			auto operating_conditions = well.control.operating_conditions.get<Flowrate_operating_conditions>();
			minimum_pressure = operating_conditions.pressure_limit;
			imposed_flowrate = operating_conditions.flowrate;
		}
		//std::cerr << "Create Producer_data_info: " << radius << " " << minimum_pressure << " " << imposed_flowrate << std::endl;
	}
};

struct Injector_data_info
{
	double radius;
	double temperature;
	double maximum_pressure;
	double imposed_flowrate;
	char operating_code;  // 'p' for pressure mode; 'f' for flowrate mode
	Injector_data_info(const Well& well) {
		assert(well.is_injecting());
		radius = well.geometry.radius;
		assert(well.control.status.is<Injection_well_status>());
		auto status = well.control.status.get<Injection_well_status>();
		temperature = status.temperature;
		if (well.operates_on_pressure()) {
			operating_code = 'p';
			assert(well.control.operating_conditions.is<Pressure_operating_conditions>());
			auto operating_conditions = well.control.operating_conditions.get<Pressure_operating_conditions>();
			maximum_pressure = operating_conditions.pressure;
			imposed_flowrate = operating_conditions.flowrate_limit;
		}
		else {
			assert(well.operates_on_flowrate());
			operating_code = 'f';
			assert(well.control.operating_conditions.is<Flowrate_operating_conditions>());
			auto operating_conditions = well.control.operating_conditions.get<Flowrate_operating_conditions>();
			maximum_pressure = operating_conditions.pressure_limit;
			imposed_flowrate = operating_conditions.flowrate;
		}
		if (imposed_flowrate > 0) {
			imposed_flowrate = -imposed_flowrate;
			std::cerr << "WARNING - setting negative injection flowrate" << std::endl;
		}
		//std::cerr << "Create Injector_data_info: " << radius << " " << temperature << " " << maximum_pressure << " " << imposed_flowrate << std::endl;
	}
};

void set_well_data(py::list wells)
{
	std::vector<Producer_data_info> producers_info;
	std::vector<Injector_data_info> injectors_info;
	for (auto&& p : wells) {
		auto well = p.cast<const Well&>();
		if (well.is_producing()) producers_info.push_back(Producer_data_info{ well });
		if (well.is_injecting()) injectors_info.push_back(Injector_data_info{ well });
	}
	producers_info.shrink_to_fit();
	injectors_info.shrink_to_fit();
	auto wrapped_producers_info = ArrayWrapper::wrap(producers_info);
	auto wrapped_injectors_info = ArrayWrapper::wrap(injectors_info);
	Well_set_wells_data(wrapped_producers_info, wrapped_injectors_info);
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
			assert(instance.control.operating_conditions.is<Pressure_operating_conditions>());
			auto operating_conditions = instance.control.operating_conditions.get<Pressure_operating_conditions>();
			return py::make_tuple(operating_conditions.pressure, operating_conditions.flowrate_limit);
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
			assert(instance.control.operating_conditions.is<Flowrate_operating_conditions>());
			auto operating_conditions = instance.control.operating_conditions.get<Flowrate_operating_conditions>();
			return py::make_tuple(operating_conditions.flowrate, operating_conditions.pressure_limit);
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

	module.def("set_well_geometries", &set_well_geometries, "Set well geometries.");
	module.def("set_well_data", &set_well_data, "Set well data.");


    py::class_<Fortran_well_data>(module, "WellData")
        .def_readwrite("operating_code", &Fortran_well_data::operating_code)
        .def_readwrite("radius", &Fortran_well_data::radius)
        .def_readwrite("maximum_pressure", &Fortran_well_data::maximum_pressure)
        .def_readwrite("minimum_pressure", &Fortran_well_data::minimum_pressure)
        .def_readwrite("imposed_flowrate", &Fortran_well_data::imposed_flowrate)
        .def_readwrite("injection_temperature", &Fortran_well_data::injection_temperature)
        ;

}
