#include <array>
#include <vector>
#include <string>

#include <pybind11/pybind11.h>
// mandatory to have optional_caster and variant_caster
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "mapbox/variant.hpp"

#include "meshtools.h"
#include "meshtools-helpers.h"

// cf. doc at http://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html?highlight=variant#c-17-library-containers
namespace pybind11 {
	namespace detail {
		template <typename... Ts>
		struct type_caster<mapbox::util::variant<Ts...>> : variant_caster<mapbox::util::variant<Ts...>> {};
		template <>
		struct visit_helper<mapbox::util::variant> {
			template <typename... Args>
			static auto call(Args &&...args) -> decltype(mapbox::util::apply_visitor(args...)) {
				return mapbox::util::apply_visitor(args...);
			}
		};
	}
} // namespace pybind11::detail

namespace py = pybind11;
namespace MT = MeshTools;

typedef MT::Point<double, 3> Point;
template <typename ... CellTypes>
using Mesh = MT::Mesh<Point, CellTypes...>;

// The following makes all STL vectors opaque
// (cf. PYBIND11_MAKE_OPAQUE)
namespace pybind11 {
	namespace detail {
		template<typename T>
		class type_caster<std::vector<T>> :
			public type_caster_base<std::vector<T>>
		{};
	}
}

template <typename TargetType, std::size_t n, typename T>
inline auto build_from_array_like(py::array_t<T, py::array::c_style > a)
-> std::unique_ptr<TargetType>
{
	typedef const std::array<T, n> Array;
	//a.attr("shape") = py::make_tuple(1);
	assert(a.ndim() == 1);
	assert(a.shape(0) == n);
	assert(a.strides(0) == a.itemsize());
	return std::make_unique<TargetType>(*(reinterpret_cast<Array*>(a.mutable_data(0))));
}

template <typename Element>
void add_element(py::module& module)
{
	py::class_<Element>(module, Element::name, py::buffer_protocol())
		.def(py::init([](py::array_t<MT::NodeId, py::array::c_style > a) {
		return build_from_array_like<Element, Element::nbnodes()>(a);
	}))
		// OPTIMIZE: As Elements have small footprint in memory that might be more efficient to perform a direct copy
		.def_buffer([](Element& self) -> py::buffer_info {
		return py::buffer_info(
			self.nodes.data(), sizeof(MT::NodeId),
			py::format_descriptor<MT::NodeId>::format(),
			1, { Element::nbnodes() }, { sizeof(MT::NodeId) }
		);
	})
		.def("__str__", [](const Element& element) {
		std::ostringstream buffer;
		buffer << Element::name << element.nodes;
		return buffer.str();
	});
}

void add_geometric_elements(py::module& module)
{
	add_element<MT::Triangle>(module);
	add_element<MT::Quad>(module);
	add_element<MT::Tetrahedron>(module);
	add_element<MT::Wedge>(module);
	add_element<MT::Hexahedron>(module);
}

template <typename Mesh>
auto add_mesh(py::module module)
{
	typedef typename Mesh::Mesh_traits Mesh_traits;
	typedef typename Mesh::Connectivity Connectivity;
	py::class_<Connectivity>(module, "Connectivity", py::module_local())
		.def_readwrite("cells", &Connectivity::cells)
		.def_readwrite("faces", &Connectivity::faces)
		.def("boundary_faces", &Connectivity::boundary_faces)
		.def("update_from_cellnodes", &Connectivity::update_from_cellnodes);
	typedef typename Connectivity::Cells Cells;
	py::bind_vector<typename Cells::Nodes>(module, "CellsNodes");
	py::bind_vector<typename Cells::Faces>(module, "CellsFaces");
	py::class_<Cells>(module, "Cells", py::module_local())
		.def_readwrite("nodes", &Cells::nodes)
		.def_readwrite("faces", &Cells::faces);
	typedef typename Connectivity::Faces Faces;
	py::bind_vector<typename Faces::Nodes>(module, "FacesNodes");
	py::class_<Faces>(module, "Faces", py::module_local())
		.def_readwrite("nodes", &Faces::nodes)
		.def_readwrite("cells", &Faces::cells)
		.def("id", [](const Faces& self, const typename Mesh_traits::Facet_type& facet) {
		return self.ids.at(Connectivity::make_face_index(facet));
	});
	py::class_<Mesh>(module, "Mesh", py::module_local())
		.def(py::init<>())
		.def_readwrite("vertices", &Mesh::vertices)
		.def_readwrite("connectivity", &Mesh::connectivity)
		.def("cell_centers", &Mesh::cell_centers)
		.def("face_centers", &Mesh::face_centers)
		.def("boundary_faces", [](const Mesh& mesh) {
		return mesh.connectivity.boundary_faces();
	})
		.def_property_readonly("nb_vertices", [](const Mesh& mesh) {
		return mesh.vertices.size();
	})
		.def_property_readonly("nb_cells", [](const Mesh& mesh) {
		return mesh.connectivity.cells.nb();
	})
		.def_property_readonly("nb_faces", [](const Mesh& mesh) {
		return mesh.connectivity.faces.nb();
	});
}

template <typename MeshType>
void set_mesh_vertices(MeshType& mesh,
	py::array_t<typename Point::Coordinate_type, py::array::c_style> vertices)
{
	typedef typename MeshType::Vertex_type Vertex;
	typedef typename MeshType::Coordinate_type Coordinate;
	constexpr auto dim = Point::dimension;
	static_assert(sizeof(Vertex) == dim * sizeof(Coordinate), "Inconsistent sizes in memory.");
	vertices.attr("shape") = py::make_tuple(-1, dim);
	auto raw = vertices.template unchecked<2>();
	mesh.vertices.clear();
	const auto n = raw.shape(0);
	mesh.vertices.reserve(n);
	auto first = reinterpret_cast<const Vertex*>(raw.data(0, 0));
	std::copy(first, first + n, std::back_inserter(mesh.vertices));
}

template <typename MeshType>
void set_uniform_mesh_cellnodes(MeshType& mesh, py::array_t<MT::NodeId, py::array::c_style> cells)
{
	typedef typename MeshType::Mesh_traits::Cell_type Cell_type;
	constexpr auto nn = Cell_type::nbnodes();
	auto& cellnodes = mesh.connectivity.cells.nodes;
	static_assert(sizeof(Cell_type) == nn * sizeof(MT::NodeId), "Inconsistent sizes in memory.");
	cells.attr("shape") = py::make_tuple(-1, nn);
	auto raw = cells.unchecked<2>();
	cellnodes.clear();
	const auto n = raw.shape(0);
	cellnodes.reserve(n);
	auto first = reinterpret_cast<const Cell_type*>(raw.data(0, 0));
	std::copy(first, first + n, std::back_inserter(cellnodes));
}

// CHECKME: output type delcaration (i.e. -> MeshType) is necessary
//          for "decltype workaround" to work in add_uniform_mesh_submodule for MSVC
//          (cf. line with  (decltype(&make_uniform_mesh<Mesh_type>))&make_uniform_mesh<Mesh_type> there)
//          "decltype workaround" itself is due to a gcc bug...
template <typename MeshType>
auto make_uniform_mesh(py::object vertices, py::object cells) -> MeshType
{
	auto mesh = MeshType{};
	set_mesh_vertices(mesh, vertices);
	set_uniform_mesh_cellnodes(mesh, cells);
	mesh.connectivity.update_from_cellnodes();
	return mesh;
}

template <typename MeshType>
auto add_mesh_submodule(py::module& module, const std::string& name)
{
	auto submodule = module.def_submodule(name.c_str());
	add_mesh<MeshType>(submodule);
	return submodule;
}

template <typename ElementType>
auto add_uniform_mesh_submodule(py::module& module, const std::string& name)
{
	typedef Mesh<ElementType> Mesh_type;
	auto submodule = add_mesh_submodule<Mesh_type>(module, name);
	// CHECKME: The weird conversion is due to a gcc bug
	// cf. https://stackoverflow.com/questions/45077622/using-a-template-function-pointer-inside-another-template-function?noredirect=1#comment77158283_45077622
	submodule.def("make", (decltype(&make_uniform_mesh<Mesh_type>))&make_uniform_mesh<Mesh_type>);
	return submodule;
}

void add_mesh_tools(py::module& module)
{

	module.doc() = "pybind11 meshtools";

	add_geometric_elements(module);

	py::class_<Point>(module, "Point")
		.def(py::init([](py::array_t<typename Point::Coordinate_type, py::array::c_style > a) {
		return build_from_array_like<Point, Point::dimension>(a);
	}))
		.def("__str__", [](const Point& P) {
		std::ostringstream buffer;
		buffer << P;
		return buffer.str();
	});

	py::class_<MT::FaceNeighbors>(module, "FaceNeighbors")
		.def("is_on_boundary", &MT::FaceNeighbors::is_on_boundary)
		.def("is_inside", &MT::FaceNeighbors::is_inside)
		.def("__str__", [](const MT::FaceNeighbors& self) {
		auto neighbors = self.value();
		std::ostringstream buffer;
		buffer << "Neighbors(" << neighbors.first << "|" << neighbors.second << ")";
		return buffer.str();
	});

	py::bind_vector<std::vector<MT::ElementId>>(module, "IDVector");
	py::bind_vector<std::vector<Point>>(module, "Vertices");
	py::bind_vector<std::vector<MT::FaceNeighbors>>(module, "FacesCells");

	add_uniform_mesh_submodule<MT::Tetrahedron>(module, "TetMesh");
	add_uniform_mesh_submodule<MT::Hexahedron>(module, "HexMesh");
	add_mesh_submodule<Mesh<MT::Tetrahedron, MT::Wedge>>(module, "TetWedgeMesh");
	add_mesh_submodule<Mesh<MT::Tetrahedron, MT::Wedge, MT::Hexahedron>>(module, "HybridMesh");

	module.def("as_coordinate_array", [](py::object object) {
		std::vector<Point>& v = object.cast<std::vector<Point>&>();
		return py::array_t<Point::Coordinate_type, py::array::c_style> {
			{v.size(), Point::dimension},
				reinterpret_cast<Point::Coordinate_type*>(v.data()->data()),
				object
		};
	}, py::keep_alive<0, 1>());
	module.def("idtype", []() { return py::dtype::of<MT::ElementId>(); });
	module.def("as_id_array", [](py::object object) {
		std::vector<MT::ElementId>& v = object.cast<std::vector<MT::ElementId>&>();
		return py::array_t<MT::ElementId, py::array::c_style> {
			{v.size()},
				reinterpret_cast<MT::ElementId*>(v.data()),
				object
		};
	}, py::keep_alive<0, 1>());

}
