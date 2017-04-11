#include <cassert>
#include <sstream>

#include "Well.h"

#include "Well_wrappers.h"

using namespace ComPASS::Well;

std::ostream& operator<<(std::ostream& os, const Well_geometry::Branch_geometry& branch)
{
	os << "Branch:(";
	std::copy(branch.begin(), branch.end(), std::ostream_iterator<Well_geometry::Node_id_type>(os, ", "));
	os << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, const Well& well)
{
	auto& geometry = well.geometry;
	os << "Well with radius " << geometry.radius << std::endl;
	os << "and " << geometry.branches.size() << " branches:";
	for (auto&& branch : geometry.branches) {
		os << std::endl << branch;
	}
	return os;
}


void add_well_wrappers(py::module& module)
{

	py::class_<Well_geometry>(module, "WellGeometry")
		.def(py::init<>())
		.def_readwrite("radius", &Well_geometry::radius)
		.def("add_branch", [](Well_geometry& instance, py::iterable& nodes) {
		Well_geometry::Branch_geometry branch;
		for (auto&& node : nodes) {
			branch.push_back(node.cast<int>());
		}
		assert(branch.size() > 1);
		instance.branches.push_back(branch);
	});

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
