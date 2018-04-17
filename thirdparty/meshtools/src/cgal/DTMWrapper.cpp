#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "DTM.h"
#include "DTMWrapper-module.h"

namespace BRGM = ExperimentalBRGM;

typedef CGAL::Epick Kernel;

typedef typename Kernel::Triangle_3 Triangle;
typedef BRGM::DTM<Triangle> DTM;
typedef DTM::Point Point;
typedef Point::FT Coordinate_type;
constexpr auto dim = Point::Ambient_dimension::value;

typedef int ShapeConnectivityId;

namespace py = pybind11;

auto dtm_from_triangles(py::array_t<Coordinate_type, py::array::c_style> vertices,
	                    py::array_t<ShapeConnectivityId, py::array::c_style> triangles)
{
	vertices.attr("shape") = py::make_tuple(-1, dim);
	triangles.attr("shape") = py::make_tuple(-1, 3);
	const auto n = triangles.shape(0);
	auto patches = DTM::Patches{};
	patches.reserve(n);
	auto tri = triangles.unchecked<2>();
	auto trivert = [&](int i, int j) { 
		return *(reinterpret_cast<const Point *>(vertices.data(tri(i, j), 0))); 
	};
	for (int i = 0; i < n; ++i) {
		patches.emplace_back(
			trivert(i, 0),
			trivert(i, 1),
			trivert(i, 2)
		);
	}
	return DTM{ patches };
}

auto triangulate_points(py::array_t<Coordinate_type, py::array::c_style> points)
{
	points.attr("shape") = py::make_tuple(-1, dim);
	auto raw = points.unchecked<2>();
	const auto n = raw.shape(0);
	auto first = reinterpret_cast<const Point *>(raw.data(0, 0));
	return BRGM::build_dtm_from_triangulation<Kernel>(first, first + n);
}

auto compute_depths(const DTM& dtm, py::array_t<Coordinate_type, py::array::c_style> points)
{
	points.attr("shape") = py::make_tuple(-1, dim);
	auto rawpts = points.unchecked<2>();
	auto n = rawpts.shape(0);
	auto P = reinterpret_cast<const Point *>(rawpts.data(0, 0));
	auto result = py::array_t<Coordinate_type, py::array::c_style>(n);
	auto rawres = result.mutable_unchecked<1>();
	auto d = reinterpret_cast<Coordinate_type *>(rawres.mutable_data(0));
	for (; n > 0; --n) {
		*d = dtm.depth(*P);
		++d;  ++P;
	}
	return result;
}

void add_dtm_wrapper(py::module& module)
{

    module.doc() = "pybind11 quick and dirty DTM wrapper";

	py::class_<DTM>(module, "DTM")
		.def("depths", &compute_depths)
		.def("nb_triangles", &DTM::number_of_shapes);

	module.def("dtm_from_triangles", &dtm_from_triangles);
	module.def("triangulate_points", &triangulate_points);

	module.def("id_dtype", []() { return py::dtype::of<ShapeConnectivityId>(); });

}
