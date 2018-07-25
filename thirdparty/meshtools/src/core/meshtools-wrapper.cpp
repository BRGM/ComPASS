#include <array>
#include <vector>
#include <algorithm>
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

namespace pyMeshTools {
    typedef int8_t byte;
    static_assert(sizeof(byte) == 1, "Inconsistent byte size.");
    typedef py::array_t<byte, py::array::c_style> RawArray;
}

namespace pyMT = pyMeshTools;

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
		.def("vtk_id", [](const Element&) {
		return Element::VTK_ID;
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
    add_element<MT::Pyramid>(module);
    add_element<MT::Hexahedron>(module);
}

template <typename Collection>
auto extract_ids(const Collection& elements) -> py::array_t<MT::VTK_ID_type, py::array::c_style>
{
	auto result = py::array_t<MT::VTK_ID_type, py::array::c_style>(elements.nb());
	auto buffer = result.request();
	auto p = reinterpret_cast<MT::VTK_ID_type*>(buffer.ptr);
	for (const auto& element : elements.nodes) {
		*p = vtk_id(element);
		++p;
	}
	return result;
}

// CHECKME: This version relies on std::map binary search to build part without knowing their number
//          when the number of part is known the other implementation might be more efficient
template <typename ItemId, typename ProcId>
auto distribute_items(const ProcId *first, const ProcId *last)
{
	std::map<ProcId, std::vector<ItemId>> item_map;
	auto id = std::size_t{ 0 };
	for (auto p = first; p != last; ++p) {
		item_map[*p].emplace_back(static_cast<ItemId>(id));
		++id;
	}
	std::vector<std::vector<ItemId>> owns;
	for (auto& part: item_map) {
		owns.emplace_back();
		owns.back().swap(part.second); // avoid copy
	}
	return owns;
}

// IMPROVE: Use hint on number of items per proc to avoid to much reallocation?
template <typename ItemId, typename ProcId>
auto distribute_items(const std::size_t nbprocs, const std::size_t nbitems, const ProcId *itemproc)
{
	auto owns = std::vector< std::vector<ItemId> >{ nbprocs };
	for (std::size_t id=0; id < nbitems; ++id) {
		assert(itemproc[id] >= 0);
		const auto proc = static_cast<std::size_t>(itemproc[id]);
		assert(proc < nbprocs);
		owns[proc].emplace_back(static_cast<ItemId>(id));
	}
	return owns;
}

struct GhostOffset {
	std::size_t cells;
	std::size_t nodes;
	std::size_t faces;
};

// FIXME: This is to replace std::vector<T>(std::size_t count, const T& value) constructor
//        which does not behave as expected with MSVC++
template <typename T>
auto filled_vector(std::size_t count, const T& value)
{
	auto result = std::vector<T>{};
	result.reserve(count);
	for (; count != 0; --count) {
		result.emplace_back(value);
	}
	return result;
}

// CHECKME: output type delcaration (i.e. -> py::list) is necessary
//          for "decltype workaround" to work for MSVC
//          cf. line with  (decltype(&distribute<Mesh>))&distribute<Mesh> there
//          "decltype workaround" itself is due to a gcc bug...
template <typename Mesh, typename ProcId=int>
auto distribute(const Mesh& mesh, py::array_t<ProcId, py::array::c_style> cellproc) -> py::list
{
	auto buffer = cellproc.request();
	assert(buffer.ndim == 1);
	assert(buffer.size == mesh.nb_cells());
	const auto first = reinterpret_cast<const ProcId *>(buffer.ptr);
	const auto last = first + buffer.size;
	auto own_cells = distribute_items<MeshTools::CellId>(first, last);
	const auto nbprocs = own_cells.size();
	// FIXME/IMPROVE: things could be factorized here if node range was acces by range and not node_range 
	// Process nodes 
	auto nodeproc = filled_vector(mesh.nb_vertices(), nbprocs);
	{
		auto pproc = first;
		for (auto&& cellnodes : mesh.connectivity.cells.nodes) {
			auto proc = *pproc;
			for (auto&& nid : MeshTools::node_range(cellnodes)) {
				if (proc < nodeproc[nid]) {
					nodeproc[nid] = proc;
				}
			}
			++pproc;
		}
		assert(pproc == last);
	}
	auto own_nodes = distribute_items<MeshTools::NodeId>(nbprocs, nodeproc.size(), nodeproc.data());
	// Process faces - very similar to process nodes (cf. improve remark)
	auto faceproc = filled_vector(mesh.nb_faces(), nbprocs);
	{
		auto pproc = first;
		for (auto&& cellfaces : mesh.connectivity.cells.faces) {
			auto proc = *pproc;
			for (auto&& fid : MeshTools::range(cellfaces)) {
				if (proc < faceproc[fid]) {
					faceproc[fid] = proc;
				}
			}
			++pproc;
		}
	}
	auto own_faces = distribute_items<MeshTools::FaceId>(nbprocs, faceproc.size(), faceproc.data());
	// End of own processing
	auto ghost_cells = std::vector< std::set<MeshTools::CellId> >{ nbprocs };
	{
		auto pproc = first;
		for (const auto& face : mesh.connectivity.faces.cells) {
			const auto pair = face.value();
			if (pproc[pair.first] != pproc[pair.second]) {
				ghost_cells[pproc[pair.first]].insert(pair.second);
				ghost_cells[pproc[pair.second]].insert(pair.first);
			}
		}
	}
	auto ghost_offsets = std::vector<GhostOffset>{};
	ghost_offsets.reserve(nbprocs);
	for (std::size_t proc = 0; proc < nbprocs; ++proc) {
		auto offset = GhostOffset{};
		auto& cells = own_cells[proc];
		auto& nodes = own_nodes[proc];
		auto& faces = own_faces[proc];
		offset.cells = cells.size();
		offset.nodes = nodes.size();
		offset.faces = faces.size();
		ghost_offsets.emplace_back(offset);
		auto& ghosts = ghost_cells[proc];
		cells.insert(cells.end(), ghosts.begin(), ghosts.end());
		// Process nodes
		{
			auto ghost_nodes = std::set<MeshTools::NodeId>{};
			for (auto&& cid : cells) {
				const auto& cellnodes = mesh.connectivity.cells.nodes[cid];
				for (auto&& nid : MeshTools::node_range(cellnodes)) {
					if (nodeproc[nid] != proc) {
						ghost_nodes.insert(nid);
					}
#ifndef NDEBUG
					else {
						assert(std::find(nodes.begin(), nodes.end(), nid) != nodes.end());
					}
#endif // NDEBUG
				}
			}
			nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
		}
		// Process faces
		{
			auto ghost_faces = std::set<MeshTools::FaceId>{};
			for (auto&& cid : cells) {
				const auto& cellfaces = mesh.connectivity.cells.faces[cid];
				for (auto&& fid : MeshTools::range(cellfaces)) {
					if (faceproc[fid] != proc) {
						ghost_faces.insert(fid);
					}
#ifndef NDEBUG
					else {
						assert(std::find(faces.begin(), faces.end(), fid) != faces.end());
					}
#endif // NDEBUG
				}
			}
			faces.insert(faces.end(), ghost_faces.begin(), ghost_faces.end());
		}
	}
	auto result = py::list{};
	for (std::size_t proc = 0; proc < nbprocs; ++proc) {
		const auto& offsets = ghost_offsets[proc];
		result.append(
			py::make_tuple(
				py::make_tuple(own_cells[proc], own_nodes[proc], own_faces[proc]),
				py::make_tuple(offsets.cells, offsets.nodes, offsets.faces)
				)
		);
	}
	return result;
}

template <typename Mesh>
auto locate_all_cell_faces(const Mesh& mesh, py::array_t<MeshTools::CellId, py::array::c_style> cells)
{
	typedef uint8_t FaceIndex;
	typedef std::pair<MeshTools::CellId, FaceIndex> Location;
	std::map<MeshTools::FaceId, Location> locations;
	assert(cells.ndim() == 1);
	auto pcell = cells.data();
	auto const end = pcell + cells.size();
	for (; pcell != end;++pcell) {		
		const auto cell = *pcell;
		auto faces = MeshTools::range(mesh.connectivity.cells.faces[cell]);
		auto k = FaceIndex{ 0 };
		for (auto&& face : faces) {
			// FIXME: we want try_emplace from C++17
			locations.insert(std::make_pair(face, Location{ cell, k }));
			assert(k < std::numeric_limits<FaceIndex>::max());
			++k;
		}
	}
	return locations;
}

// CHECKME: output type delcaration (i.e. -> py::tuple) is necessary
//          for "decltype workaround" to work for MSVC
//          cf. line with  (decltype(&locate_cell_faces<Mesh>))&locate_cell_faces<Mesh> there
//          "decltype workaround" itself is due to a gcc bug...
//FIXME: g++ complains with FaceIndex as a template because of
//		 "auto face_position = position.unchecked<1>();"
//       we impose FaceIndex = uint8_t (cf. below)
template <typename Mesh, typename FaceIndex = uint8_t>
auto identify_faces_from_positions(
	const Mesh& mesh,
	py::array_t<MeshTools::CellId, py::array::c_style> cell,
//	py::array_t<FaceIndex, py::array::c_style> position) // FIXME: cf.above
	py::array_t<uint8_t, py::array::c_style> position)
	-> py::array_t<MT::FaceId, py::array::c_style>
{
	assert(cell.ndim() == 1);
	assert(position.ndim() == 1);
	assert(cell.shape(0) == position.shape(0));
	const auto nf = position.shape(0);
	auto result = py::array_t<MT::FaceId, py::array::c_style>(nf);
	auto face_id = result.mutable_unchecked<1>();
	auto face_cell = cell.unchecked<1>();
	auto face_position = position.unchecked<1>();
	using Accessor = MT::VariantAccessor<MT::ArrayBaseAccessor>;
	for (std::size_t k = 0; k != nf; ++k) {
		face_id(k) = Accessor::extract(mesh.connectivity.cells.faces[face_cell(k)], face_position(k));
	}
	return result;
}

// CHECKME: output type delcaration (i.e. -> py::tuple) is necessary
//          for "decltype workaround" to work for MSVC
//          cf. line with  (decltype(&locate_cell_faces<Mesh>))&locate_cell_faces<Mesh> there
//          "decltype workaround" itself is due to a gcc bug...
template <typename Mesh>
auto locate_faces_with_cell(const Mesh& mesh, py::array_t<MeshTools::CellId, py::array::c_style> cells, py::array_t<MeshTools::FaceId, py::array::c_style> faces) -> py::tuple
{
	auto locations = locate_all_cell_faces(mesh, cells);
	typedef decltype(locations) LocationMap;
	typedef typename LocationMap::mapped_type Location;
	typedef decltype(Location::second) FaceIndex;
	assert(faces.ndim() == 1);
	const std::size_t nbfaces = faces.size();
	auto face_cell = py::array_t<MeshTools::CellId, py::array::c_style>{ nbfaces };
	auto p_face_cell = face_cell.mutable_data();
	auto face_index = py::array_t<FaceIndex, py::array::c_style>{ nbfaces };
	auto p_face_index = face_index.mutable_data();
	auto p_face = faces.data();
	const auto end = p_face + nbfaces;
	for (; p_face != end; ++p_face) {
		const auto location = locations.at(*p_face);
		(*p_face_cell) = location.first;
		++p_face_cell;
		(*p_face_index) = location.second;
		++p_face_index;
	}
	return py::make_tuple(face_cell, face_index);
}

template <typename Vector>
struct CheckBinds
{
private:
    CheckBinds() = delete;
    CheckBinds(const CheckBinds&) = delete;
    static bool already_binded;
public:
    static bool needs_binding() {
        if (already_binded) return false;
        already_binded = true;
        return true;
    }
};

template <typename Vector>
bool CheckBinds<Vector>::already_binded = false;

// CHECKME: output type declaration (i.e. -> Vector) is necessary in MSVC
//          for "decltype workaround" to work below
//          "decltype workaround" itself is due to a gcc bug...
template <typename Vector>
auto from_raw_array(pyMT::RawArray raw_array) -> Vector
{
    typedef typename Vector::value_type value_type;
    assert(raw_array.ndim() == 2);
    assert(raw_array.shape(1) == sizeof(value_type));
    auto begin = reinterpret_cast<const value_type *>(raw_array.data());
    auto end = begin + raw_array.shape(0);
    return Vector{ begin, end };
}

// CHECKME: output type declaration (i.e. -> py::array_t<...>) is necessary in MSVC
//          for "decltype workaround" to work below
//          "decltype workaround" itself is due to a gcc bug...
template <typename Vector>
auto as_raw_array(const Vector& v) -> pyMT::RawArray
{
    typedef typename Vector::value_type value_type;
    typedef pyMT::byte byte;
    return pyMT::RawArray{
        { v.size(), sizeof(value_type) }, // shape
        { sizeof(value_type), sizeof(byte) }, // stride
            reinterpret_cast<const byte*>(v.data())
    };
}

template <typename Vector>
auto pybind_vector(py::module module, const char * classname)
{
    typedef typename Vector::value_type value_type;
    typedef int8_t byte;
    static_assert(sizeof(byte) == 1, "Inconsistent byte size.");
    return py::bind_vector<Vector>(module, classname)
        .def_static("from_raw_array", (decltype(&from_raw_array<Vector>)) &from_raw_array<Vector>)  // decltype is due to gcc bug
        .def("raw_array", (decltype(&as_raw_array<Vector>)) &as_raw_array<Vector>,  // decltype is due to gcc bug
            py::keep_alive<0, 1>())
        ;
}

template <typename Vector>
void bind_vector(py::module module, const char * classname)
{
    if(CheckBinds<Vector>::needs_binding()) {
        pybind_vector<Vector>(module, classname);
    }
}

template <typename Vector>
void bind_vector_with_array_view(py::module module, const char * classname)
{
    if(CheckBinds<Vector>::needs_binding()) {
        typedef typename Vector::value_type value_type;
        pybind_vector<Vector>(module, classname).def("array_view", [](const Vector& self) {
            return py::array_t<value_type, py::array::c_style> {
                {self.size()}, { sizeof(value_type) }, self.data()
            };
        }, py::keep_alive<0, 1>());
    }
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
	typedef typename Cells::Cell_type Cell_type;
	// The following would give acess to buffer view of CellNodes vector but desactivate variant caster
	// py::class_<typename Cells::Cell_type>(module, "Cell_type", py::module_local());
	bind_vector<typename Cells::Nodes>(module, "CellsNodes");
	bind_vector<typename Cells::Faces>(module, "CellsFaces");
	py::class_<Cells>(module, "Cells", py::module_local())
		.def_readwrite("nodes", &Cells::nodes)
		.def_readwrite("faces", &Cells::faces);
	typedef typename Connectivity::Faces Faces;
    // CHECKME: Without the CheckBinds counter the following bind_vector would fail at import with 
    //          ImportError: generic_type: type "FacesNodes" is already registered!
    bind_vector<typename Faces::Nodes>(module, "FacesNodes");
    py::class_<Faces>(module, "Faces", py::module_local())
        .def_readwrite("nodes", &Faces::nodes)
        .def_readwrite("cells", &Faces::cells)
		.def("cells_as_array", [](const Faces& self) {
		return py::array_t<MeshTools::CellId, py::array::c_style>{
			std::vector<std::size_t>{ { self.cells.size(), 2 } }, // shape 
			std::vector<std::size_t>{ { 2 * sizeof(MeshTools::CellId), sizeof(MeshTools::CellId) } }, // strides 
			reinterpret_cast<const MeshTools::CellId *>(self.cells.data())
		};
	}, py::keep_alive<0, 1>() )
		.def("id", [](const Faces& self, const typename Mesh_traits::Facet_type& facet) {
		return self.ids.at(Connectivity::make_face_index(facet));
	});
    auto mesh_class =
	py::class_<Mesh>(module, "Mesh", py::module_local())
		.def(py::init<>())
		.def_readwrite("vertices", &Mesh::vertices)
		.def_readwrite("connectivity", &Mesh::connectivity)
		.def("cell_centers", &Mesh::cell_centers)
		.def("face_centers", &Mesh::face_centers)
		.def("boundary_cells", [](const Mesh& self) {
		return self.connectivity.boundary_cells();
	})
		.def("boundary_faces", [](const Mesh& self) {
		return self.connectivity.boundary_faces();
	})
		.def("vertices_array", [](const Mesh& self) {
		typedef typename Mesh::Coordinate_type Coordinate_type;
		return py::array_t<Coordinate_type, py::array::c_style>{
				std::vector<std::size_t>{ { self.vertices.size(), Mesh::dimension } }, // shape 
				std::vector<std::size_t>{ { Mesh::dimension * sizeof(Coordinate_type), sizeof(Coordinate_type) } }, // strides 
				reinterpret_cast<const Coordinate_type *>(self.vertices.data())
		};
	}, py::keep_alive<0, 1>())
		.def_property_readonly("nb_vertices", [](const Mesh& mesh) {
		return mesh.vertices.size();
	})
		.def_property_readonly("nb_cells", [](const Mesh& mesh) {
		return mesh.connectivity.cells.nb();
	})
		.def_property_readonly("nb_faces", [](const Mesh& mesh) {
		return mesh.connectivity.faces.nb();
	})
		.def("cells_nodes_as_COC", [](const Mesh& mesh) {
		return MT::node_collection_as_COC_data(mesh.connectivity.cells.nodes);
	})
		.def("faces_nodes_as_COC", [](const Mesh& mesh) {
		return MT::node_collection_as_COC_data(mesh.connectivity.faces.nodes);
	})
		.def("cells_faces_as_COC", [](const Mesh& mesh) {
		return MT::face_collection_as_COC_data(mesh.connectivity.cells.faces);
	})
		.def("COC_data", [](const Mesh& mesh) {
		return py::make_tuple(
			MT::node_collection_as_COC_data(mesh.connectivity.cells.nodes),
			MT::face_collection_as_COC_data(mesh.connectivity.cells.faces),
			MT::node_collection_as_COC_data(mesh.connectivity.faces.nodes)
		);
	})
		.def("cells_vtk_ids", [](const Mesh& mesh) {
		return extract_ids(mesh.connectivity.cells);
	})
		.def("faces_vtk_ids", [](const Mesh& mesh) {
		return extract_ids(mesh.connectivity.faces);
	})
		.def("distribute", (decltype(&distribute<Mesh>))&distribute<Mesh>) // decltype is due to gcc bug
		.def("locate_faces_with_cell", (decltype(&locate_faces_with_cell<Mesh>))&locate_faces_with_cell<Mesh>)  // decltype is due to gcc bug
		.def("identify_faces_from_positions", (decltype(&identify_faces_from_positions<Mesh>))&identify_faces_from_positions<Mesh>)  // decltype is due to gcc bug
        .def("set_vertices", &set_mesh_vertices<Mesh>)
        .def(py::pickle(
            [](const Mesh &mesh) { // __getstate__
                return py::make_tuple(
                    as_raw_array(mesh.vertices),
                    as_raw_array(mesh.connectivity.cells.nodes)
                );
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid pickled state!");
                auto mesh = Mesh{};
                mesh.vertices = from_raw_array<typename Mesh::Vertices>(t[0].cast<pyMT::RawArray>());
                mesh.connectivity.cells.nodes = from_raw_array<typename Cells::Nodes>(t[1].cast<pyMT::RawArray>());
                mesh.connectivity.update_from_cellnodes();
                return mesh;
            } ))
    ;

    return mesh_class;

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

template <typename MeshType, std::size_t N, bool ok>
struct CellTypeAttempt
{
	typedef typename MeshType::Mesh_traits::Cell_types Cell_types;
	typedef typename MeshType::Mesh_traits::Cell_type Cell_type;
	static Cell_type cast(py::object object)
	{
		std::ostringstream buffer;
		buffer << "could not convert ";
		buffer << py::str(object.attr("__class__"));
		buffer << " to mesh element";
		throw py::cast_error(buffer.str());
		assert(false); // never reached
		return Cell_type{};
	}
};

template <typename MeshType, std::size_t N>
struct CellTypeAttempt<MeshType, N, true>
{
	typedef typename MeshType::Mesh_traits::Cell_types Cell_types;
	typedef typename MeshType::Mesh_traits::Cell_type Cell_type;
	static Cell_type cast(py::object object)
	{
		typedef typename std::tuple_element<N, Cell_types>::type Attempt;
		try {
			auto cell = object.cast<Attempt>();
			return Cell_type{ cell };
		}
		catch (py::cast_error)
		{
			return CellTypeAttempt<MeshType, N+1, N+1<std::tuple_size<Cell_types>::value>::cast(object);
		}
		assert(false); // never reached
		return Cell_type{};
	}
};

template <typename MeshType>
inline auto cast_as_celltype(py::object object) -> typename MeshType::Mesh_traits::Cell_type
{
	typedef typename MeshType::Mesh_traits::Cell_types Cell_types;
	static_assert(std::tuple_size<Cell_types>::value > 0, "Empty tuple!");
	return CellTypeAttempt<MeshType, 0, 0<std::tuple_size<Cell_types>::value>::cast(object);
}

template <typename MeshType>
void set_generic_mesh_cellnodes(MeshType& mesh, py::list cells)
{
	auto& cellnodes = mesh.connectivity.cells.nodes;
	cellnodes.clear();
	const auto n = py::len(cells);
	cellnodes.reserve(n);
	for (std::size_t i = 0; i!=n;++i) {
		cellnodes.push_back(cast_as_celltype<MeshType>(cells[i]));
	}
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

// CHECKME: output type delcaration (i.e. -> MeshType) is necessary
//          for "decltype workaround" to work in add_uniform_mesh_submodule for MSVC
//          (cf. line with  (decltype(&make_uniform_mesh<Mesh_type>))&make_uniform_mesh<Mesh_type> there)
//          "decltype workaround" itself is due to a gcc bug...
template <typename MeshType>
auto make_generic_mesh(py::object vertices, py::object cells) -> MeshType
{
	auto mesh = MeshType{};
	set_mesh_vertices(mesh, vertices);
	set_generic_mesh_cellnodes(mesh, cells);
	mesh.connectivity.update_from_cellnodes();
	return mesh;
}

// CHECKME: output type delcaration (i.e. -> MeshType) is necessary
//          for "decltype workaround" to work in add_uniform_mesh_submodule for MSVC
//          (cf. line with  (decltype(&rebuild_mesh<Mesh_type>))&rebuild_mesh<Mesh_type> there)
//          "decltype workaround" itself is due to a gcc bug...
template <typename MeshType, typename byte=int8_t>
auto create_from_remap(
	py::object vertices,
	py::array_t<byte, py::array::c_style> raw_cellnodes,
	py::array_t<MT::NodeId, py::array::c_style> global_node_ids) -> MeshType
{
	assert(global_node_ids.ndim() == 1);
	const auto nn = global_node_ids.shape(0);
	auto ni = MT::NodeId{ 0 };
	auto gni = global_node_ids.unchecked<1>();
	auto new_node_ids = std::map<MT::NodeId, MT::NodeId>{};
	for (std::size_t k = 0; k != nn; ++k) {
		new_node_ids[gni(k)] = ni;
		++ni;
	}
	typedef typename MeshType::Connectivity::Cell_type Cell_type;
	assert(raw_cellnodes.ndim() == 2);
	assert(raw_cellnodes.shape(1) == sizeof(Cell_type));
	auto first = reinterpret_cast<const Cell_type*>(raw_cellnodes.data());
	auto last = first + raw_cellnodes.shape(0);
	auto mesh = MeshType{};
	set_mesh_vertices(mesh, vertices);
	assert(mesh.vertices.size() == nn);
	mesh.connectivity.remap_cellnodes(first, last, new_node_ids);
	return mesh;
}
	
template <typename MeshType>
auto add_mesh_submodule(py::module& module, const std::string& name)
{
	auto submodule = module.def_submodule(name.c_str());
	auto mesh_class = add_mesh<MeshType>(submodule);
	submodule.def("create", (decltype(&make_generic_mesh<MeshType>))&make_generic_mesh<MeshType>);
	submodule.def("create_from_remap", (decltype(&create_from_remap<MeshType>))&create_from_remap<MeshType>);
	return std::make_tuple(submodule, mesh_class);
}

template <typename ElementType>
auto add_uniform_mesh_submodule(py::module& module, const std::string& name)
{
	typedef Mesh<ElementType> Mesh_type;
    // FIXME: TO be replaced by structured bindings in C++17 (e.g: auto [submodule, mesh_class] = add_mesh_submodule...)
	auto tuple = add_mesh_submodule<Mesh_type>(module, name);
    auto submodule = std::get<0>(tuple);
    auto mesh_class = std::get<1>(tuple);
    // CHECKME: The weird conversion is due to a gcc bug
	// cf. https://stackoverflow.com/questions/45077622/using-a-template-function-pointer-inside-another-template-function?noredirect=1#comment77158283_45077622
    submodule.def("make", (decltype(&make_uniform_mesh<Mesh_type>))&make_uniform_mesh<Mesh_type>);
    mesh_class.def("face_id", [](const Mesh_type& self, const typename Mesh_type::Facet_type& facet) {
        const auto& faces_id_map = self.connectivity.faces.ids;
        auto position = faces_id_map.find(self.connectivity.make_face_index(facet));
        assert(self.nb_faces() < std::numeric_limits<MT::FaceId>::max());
        if (position == faces_id_map.end()) return static_cast<MT::FaceId>(self.nb_faces());
        return position->second;
    });
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
		.def("as_tuple", [](const MT::FaceNeighbors& self) {
		return py::make_tuple(self.value().first, self.value().second);
	})
		.def("__str__", [](const MT::FaceNeighbors& self) {
		auto neighbors = self.value();
		std::ostringstream buffer;
		buffer << "Neighbors(" << neighbors.first << "|" << neighbors.second << ")";
		return buffer.str();
	});

	//bind_vector_with_array_view<std::vector<std::size_t>>(module, "COCPointersVector");
	bind_vector_with_array_view<std::vector<MT::ElementId>>(module, "IDVector");
    bind_vector<std::vector<Point>>(module, "Vertices");
    bind_vector<std::vector<MT::Quad>>(module, "QuadVector");
    bind_vector<std::vector<MT::FaceNeighbors>>(module, "FacesCells");

    add_uniform_mesh_submodule<MT::Triangle>(module, "TSurf");
    add_uniform_mesh_submodule<MT::Tetrahedron>(module, "TetMesh");
    add_uniform_mesh_submodule<MT::Hexahedron>(module, "HexMesh");
    //add_uniform_mesh_submodule<MT::Pyramid>(module, "PyramidMesh");
    add_mesh_submodule<Mesh<MT::Tetrahedron, MT::Wedge>>(module, "TetWedgeMesh");
	// FIXME: The following would require to define 2D elements facets (as 3D objects or 1D ?)
	//add_mesh_submodule<Mesh<MT::Triangle, MT::Quad, MT::Tetrahedron, MT::Wedge, MT::Hexahedron>>(module, "HybridDimensionMesh");
    //add_mesh_submodule<Mesh<MT::Tetrahedron, MT::Wedge, MT::Hexahedron>>(module, "HybridMesh");
    add_mesh_submodule<Mesh<MT::Tetrahedron, MT::Wedge, MT::Pyramid, MT::Hexahedron>>(module, "HybridMesh");

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
