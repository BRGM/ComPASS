#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "meshtools.h"

namespace py = pybind11;
namespace MT = MeshTools;

typedef MT::ElementId ElementId;
typedef py::array_t<ElementId, py::array::c_style> IdArray;

template <typename Mesh>
void set_mesh_vertices(Mesh& mesh,
	py::array_t<typename Mesh::Coordinate, py::array::c_style> vertices)
{
	typedef typename Mesh::Vertex Vertex;
	typedef typename Mesh::Coordinate Coordinate;
	constexpr auto dim = std::tuple_size<Vertex>::value;
	static_assert(sizeof(Vertex) == dim * sizeof(Coordinate), "Inconsistent sizes in memory.");
	vertices.attr("shape") = py::make_tuple(-1, dim);
	auto raw = vertices.template unchecked<2>();
	mesh.vertices.clear();
	const auto n = raw.shape(0);
	mesh.vertices.reserve(n);
	auto first = reinterpret_cast<const Vertex*>(raw.data(0, 0));
	std::copy(first, first + n, std::back_inserter(mesh.vertices));
}

template <typename Mesh>
void set_mesh_cellnodes(Mesh& mesh,
	IdArray cells)
{
	typedef typename Mesh::Element Element;
	constexpr auto nn = Element::nb_nodes();
	auto& cellnodes = mesh.connectivity.cells.nodes;
	typedef typename std::remove_reference<decltype(cellnodes)>::type::value_type CellNodes;
	static_assert(std::tuple_size<CellNodes>::value == nn, "Inconsistent sizes.");
	static_assert(sizeof(CellNodes) == nn * sizeof(ElementId), "Inconsistent sizes in memory.");
	cells.attr("shape") = py::make_tuple(-1, nn);
	auto raw = cells.unchecked<2>();
	cellnodes.clear();
	const auto n = raw.shape(0);
	cellnodes.reserve(n);
	auto first = reinterpret_cast<const CellNodes*>(raw.data(0, 0));
	std::copy(first, first + n, std::back_inserter(cellnodes));
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto make_mesh(py::array_t<typename Mesh::Coordinate, py::array::c_style> vertices,
	IdArray cells) -> Mesh
{
	auto mesh = Mesh{};
	set_mesh_vertices(mesh, vertices);
	set_mesh_cellnodes(mesh, cells);
	mesh.connectivity.update_from_cellnodes();
	return mesh;
}

template <typename Mesh>
auto prepare_point_result(std::size_t n)
{
	typedef typename Mesh::Vertex Vertex;
	typedef typename Mesh::Coordinate Coordinate;
	constexpr auto dim = std::tuple_size<Vertex>::value;
	static_assert(sizeof(Coordinate)*dim == sizeof(Vertex), "Inconsistent sizes in memory.");
	auto result = py::array_t<Coordinate>{ n * dim };
	result.attr("shape") = py::make_tuple(n, dim);
	return result;
}

/** IMPROVE: add, strides, other or use library... and Iterator concepts validity */
template <typename Position = std::size_t>
class Range {
protected:
	struct iterator {
		Position pos;
		bool operator<(const iterator& other) const {
			return pos < other.pos;
		}
		bool operator!=(const iterator& other) const {
			return pos != other.pos;
		}
		bool operator==(const iterator& other) const {
			return pos == other.pos;
		}
		iterator& operator++() {
			++pos;
			return *this;
		}
		Position operator*() {
			return pos;
		}
		iterator(Position p) :
			pos(p) {}
	};
	Position first;
	Position last;
public:
	Range(Position p) :
		first{ 0 },
		last{ p }
	{}
	Range(Position p1, Position p2) :
		first{ p1 },
		last{ p2 }
	{}
	auto begin() const { return iterator(first); }
	auto end() const { return iterator(last); }
	auto size() const { return std::max(0, last - first); }
};

template <typename Mesh, typename IdIterator, typename ConnectivityTable>
auto compute_centers(const Mesh& mesh,
	std::size_t n, IdIterator pid,
	const ConnectivityTable& nodes)
{
	auto result = prepare_point_result<Mesh>(n);
	auto rawres = result.template mutable_unchecked<2>();
	typedef typename Mesh::Vertex Vertex;
	auto p = reinterpret_cast<Vertex*>(rawres.mutable_data(0, 0));
	for (; n != 0; --n) {
		*p = MT::compute_center(MT::extract(mesh.vertices, nodes[*pid]));
		++pid;
		++p;
	}
	return result;
}

template <typename Mesh, typename ConnectivityTable>
auto compute_centers(const Mesh& mesh,
	IdArray ids,
	const ConnectivityTable& nodes)
{
	auto rawids = ids.unchecked<1>();
	return compute_centers(mesh, rawids.size(), rawids.data(0), nodes);
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto face_centers(const Mesh& mesh, IdArray faces) -> py::array_t<typename Mesh::Coordinate>
{
	return compute_centers(mesh, faces, mesh.connectivity.faces.nodes);
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto cell_centers(const Mesh& mesh, IdArray cells) -> py::array_t<typename Mesh::Coordinate>
{
	return compute_centers(mesh, cells, mesh.connectivity.cells.nodes);
}

template <typename Mesh, typename GeometricFeature>
auto all_centers(const Mesh& mesh, const GeometricFeature& feature)
{
	const auto n = feature.nb();
	return compute_centers(mesh, n, Range<>{ n }.begin(), feature.nodes);
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto all_face_centers(const Mesh& mesh) -> py::array_t<typename Mesh::Coordinate>
{
	return all_centers(mesh, mesh.connectivity.faces);
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto all_cell_centers(const Mesh& mesh) -> py::array_t<typename Mesh::Coordinate>
{
	return all_centers(mesh, mesh.connectivity.cells);
}

// FIXME: Returned type is mandatory for MSVC and decltype conversion in add_mesh
template <typename Mesh>
auto faces_ids(const Mesh& mesh, IdArray nodes) ->  py::array_t<ElementId>
{
	typedef typename Mesh::Element Element;
	constexpr auto nfn = Element::nb_facet_nodes();
	nodes.attr("shape") = py::make_tuple(-1, nfn);
	auto rawnodes = nodes.unchecked<2>();
	auto n = rawnodes.shape(0);
	auto result = py::array_t<ElementId>{ n };
	auto rawres = result.mutable_unchecked<1>();
	typedef typename Mesh::Connectivity::Faces::Index FIndex;
	typedef typename FIndex::NodesList FaceNodes;
	static_assert(std::tuple_size<FaceNodes>::value == nfn, "Inconsistent sizes.");
	auto p = rawres.mutable_data(0);
	auto fnodes = reinterpret_cast<const FaceNodes*>(rawnodes.data(0, 0));
	const auto& facesids = mesh.connectivity.faces.ids;
	for (; n != 0; --n) {
		(*p) = facesids.at(FIndex{ *fnodes });
		++p;
		++fnodes;
	}
	return result;
}

template <typename Mesh>
auto boundary_faces(const Mesh& mesh) -> py::array_t<ElementId>
{
	// FIXME: This copy is superfluous (use py::array_t and keep_alive policy)
	auto boundaries = mesh.connectivity.boundary_faces();
	const auto n = boundaries.size();
	auto result = py::array_t<ElementId>{ n };
	auto rawres = result.mutable_unchecked<1>();
	std::copy(boundaries.begin(), boundaries.end(), rawres.mutable_data(0));
	return result;
}

template <typename Mesh>
void add_mesh(py::module& module, const char *classname, const char *factoryname)
{
	// FIXME: The decltype conversions are necessary for gcc
	// cf. https://stackoverflow.com/questions/45077622/using-a-template-function-pointer-inside-another-template-function
	// FIXME/IMPROVE: use overload cast!
	py::class_<Mesh>(module, classname)
		.def("nb_vertices", [](const Mesh& mesh) { return mesh.vertices.size(); })
		.def("nb_cells", [](const Mesh& mesh) { return mesh.connectivity.cells.nb(); })
		.def("nb_faces", [](const Mesh& mesh) { return mesh.connectivity.faces.nb(); })
		.def_property_readonly("vertices", [](const Mesh& mesh) {
		constexpr auto dim = std::tuple_size<typename Mesh::Vertex>::value;
		const auto& vertices = mesh.vertices;
		return py::array_t<typename Mesh::Coordinate, py::array::c_style>{
			{ vertices.size(), dim }, { dim * sizeof(typename Mesh::Coordinate), sizeof(typename Mesh::Coordinate) }, vertices.data()->data()
		};
	})
		.def_property_readonly("cellfaces", [](const Mesh& mesh) {
		constexpr auto nf = Mesh::Element::nb_facets();
		const auto& cells = mesh.connectivity.cells;
		return py::array_t<ElementId, py::array::c_style>{
			{ cells.nb(), nf }, { nf * sizeof(ElementId), sizeof(ElementId) }, cells.faces.data()->data()
		};
	})
		.def_property_readonly("cellnodes", [](const Mesh& mesh) {
		constexpr auto nn = Mesh::Element::nb_nodes();
		const auto& cells = mesh.connectivity.cells;
		return py::array_t<ElementId, py::array::c_style>{
			{ cells.nb(), nn }, { nn * sizeof(ElementId), sizeof(ElementId) }, cells.nodes.data()->data()
		};
	})
		.def_property_readonly("facenodes", [](const Mesh& mesh) {
		constexpr auto nfn = Mesh::Element::nb_facet_nodes();
		const auto& faces = mesh.connectivity.faces;
		return py::array_t<ElementId, py::array::c_style>{
			{ faces.nb(), nfn }, { nfn * sizeof(ElementId), sizeof(ElementId) }, faces.nodes.data()->data()
		};
	})
		//.def("cell_centers", py::overload_cast<const Mesh&, IdArray>(&cell_centers<Mesh>))
		//.def("cell_centers", py::overload_cast<const Mesh&>(&all_cell_centers<Mesh>))
		//.def("face_centers", py::overload_cast<const Mesh&, IdArray>(&face_centers<Mesh>))
		//.def("face_centers", py::overload_cast<const Mesh&>(&all_face_centers<Mesh>))
		.def("cell_centers", (decltype(&cell_centers<Mesh>))&cell_centers<Mesh>)
		.def("all_cell_centers", (decltype(&all_cell_centers<Mesh>))&all_cell_centers<Mesh>)
		.def("face_centers", (decltype(&face_centers<Mesh>))&face_centers<Mesh>)
		.def("all_face_centers", (decltype(&all_face_centers<Mesh>))&all_face_centers<Mesh>)
		.def("faces_ids", (decltype(&faces_ids<Mesh>))&faces_ids<Mesh>)
		.def("boundary_faces", (decltype(&boundary_faces<Mesh>))&boundary_faces<Mesh>);
	module.def(factoryname, (decltype(&make_mesh<Mesh>))&make_mesh<Mesh>);
}

py::module& add_mesh_tools(py::module& module)
{
	add_mesh<MT::TetMesh>(module, "TetMesh", "tetmesh");
	module.def("idtype", []() { return py::dtype::of<ElementId>(); });
	return module;
}
