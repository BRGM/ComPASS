#pragma once

#include <cassert>
#include <array>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <utility>
#include <iterator>

namespace MeshTools
{

	/** Elements Connectivity is exactly the same as in the VTK library.
	cf. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf */

	typedef int ElementId;
	typedef ElementId NodeId;
	typedef ElementId CellId;
	typedef ElementId FaceId;
	typedef int VTK_ID_type;
	template <typename FacetElementType, ::std::size_t N>
	using FacetList = ::std::array< const ::std::array<const NodeId, FacetElementType::nb_nodes>, N>;

	struct Vertex_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 1 };
		static constexpr auto facets = ::std::array< ::std::array<NodeId, 0>, 0>{};
		//static constexpr FacetList<DummyElement, 0>
		//	facets{};
		static constexpr auto nb_facet_nodes = ::std::size_t{ 0 };
		static constexpr auto VTK_ID = VTK_ID_type{ 1 };
	};

	struct Segment_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 2 };
		static constexpr auto facets = FacetList<Vertex_info, 2>{ { { 0 }, { 1 } } };
		static constexpr auto nb_facet_nodes = ::std::size_t{ 1 };
		static constexpr auto VTK_ID = VTK_ID_type{ 3 };
	};

	struct Triangle_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 3 };
		// Facets are CGAL compliant : Face n correspond to nth vertex missing
		static constexpr auto facets = FacetList<Segment_info, 3>{ { { 1, 2 },{ 0, 2 },{ 0, 1 } } };
		static constexpr auto nb_facet_nodes = ::std::size_t{ 2 };
		static constexpr auto VTK_ID = VTK_ID_type{ 5 };
	};

	struct Quad_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 4 };
		static constexpr auto facets = FacetList<Segment_info, 4>{ { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 0, 3 } } };
		static constexpr auto nb_facet_nodes = ::std::size_t{ 2 };
		static constexpr auto VTK_ID = VTK_ID_type{ 9 };
	};

	struct Tetrahedron_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 4 };
		// CHECKME: outward oriented
		// Facets are CGAL compliant : Face n correspond to nth vertex missing
		static constexpr auto facets = FacetList<Triangle_info, 4>{ { { 1, 2, 3 },{ 0, 3, 2 },{ 0, 1, 3 },{ 0, 2, 1 } } };
		static constexpr auto nb_facet_nodes = ::std::size_t{ 3 };
		static constexpr auto VTK_ID = VTK_ID_type{ 10 };
	};

	struct Hexahedron_info {
		// CHECKME: to be removed cf. compiler error below
		static constexpr auto nb_nodes = ::std::size_t{ 8 };
		// CHECKME: outward oriented
		static constexpr auto facets = FacetList<Quad_info, 6>{ { { 0, 1, 2, 3 },
		                                                          { 4, 5, 6, 7 },
																  { 1, 2, 6, 5 },
																  { 2, 6, 7, 3 },
																  { 3, 7, 4, 0 },
																  { 0, 1, 5, 4 } } };
		static constexpr auto nb_facet_nodes = ::std::size_t{ 4 };
		static constexpr auto VTK_ID = VTK_ID_type{ 12 };
	};

	template <typename ElementType>
	constexpr auto max_node_id()
	{
		typedef typename decltype(ElementType::facets)::value_type facet_type;
		return ::std::accumulate(
			ElementType::facets.begin(), ElementType::facets.end(), 0,
			[](NodeId result, const facet_type& facet) {
			return ::std::max(result, facet.empty() ? 0 : *::std::max_element(facet.begin(), facet.end()));
		}
		);
	}

	//struct Node {};
	//struct Edge {};
	//struct Face {};
	//struct Cell {};

	template <typename ElementInfo>
	class Element {
		//public:
	protected:
		// CHECKME: Compiler error???? The following does not work
		// static constexpr auto nn = max_node_id<ElementInfo>() + 1;
		static constexpr auto nn = ElementInfo::nb_nodes;
		static constexpr auto nf = ElementInfo::facets.size();
		static constexpr auto nfn = ElementInfo::nb_facet_nodes; // FIXME: this supposes identical facets
	public:
		typedef ElementInfo info;
		static constexpr auto nb_nodes() { return nn; }
		static constexpr auto nb_facets() { return nf; }
		static constexpr auto nb_facet_nodes() { return nfn; }
		static constexpr auto vtk_id() { return ElementInfo::VTK_ID; }
		//template <typename Nodes>
		//static auto project(const Nodes& nodes) {
		//	::std::array<Nodes::value_type, nf>
		//}
		template<typename CellNodes, typename FacetNodes>
		static void facet_nodes(const CellNodes& cellnodes, int facet, FacetNodes& facet_nodes) {
			assert(facet >= 0 && facet < nf);
			const auto& local_nodes = info::facets[facet];
			for (::std::size_t i = 0; i < nfn; ++i) {
				facet_nodes[i] = cellnodes[local_nodes[i]];
			}
		}
		template<typename CellNodes>
		static auto facet_nodes(const CellNodes& cellnodes, int facet) {
			assert(facet >= 0 && facet < nf);
			auto facet_nodes = info::facets[facet];
			facet_nodes(cellnodes, facet, facet_nodes);
			return facet_nodes;
		}
		template<typename CellNodes, typename OutputIterator>
		static auto collect_facet_nodes(const CellNodes& cellnodes, OutputIterator out) {
			for (int facet = 0; facet < nf; ++facet) {
				facet_nodes(cellnodes, facet, *(++out));
			}
			return out;
		}
		//template <typename LocusType>
		//struct Nb {};
		//template<>
		//struct Element<ElementInfo>::Nb<Node> {
		//	static constexpr auto value = Element<ElementInfo>::nb_nodes();
		//};
		//template<>
		//struct Element<ElementInfo>::Nb<Face> {
		//	static constexpr auto value = Element<ElementInfo>::nb_facets();
		//};

	};


	typedef Element<Vertex_info> Vertex;
	typedef Element<Segment_info> Segment;
	typedef Element<Triangle_info> Triangle;
	typedef Element<Quad_info> Quad;
	typedef Element<Tetrahedron_info> Tetrahedron;
	typedef Element<Hexahedron_info> Hexahedron;

	using Tet = Tetrahedron;
	using Hex = Hexahedron;

	template <typename T, ::std::size_t N>
	using FSCoC = ::std::vector< ::std::array<T, N> >;

	template <typename RowType, typename ColumnType, ::std::size_t N>
	struct HomogeneousConnectivityTable :
		FSCoC<NodeId, N> {
		typedef RowType row_type;
		typedef ColumnType col_type;
	};

	//template <typename ElementType>
	//struct LightweightHomogeneousMesh {
	//	typedef ElementType element_type;
	//	struct Cells {
	//		HomogeneousConnectivityTable<Cell, Node, ElementType::nb_nodes()> nodes;
	//	};
	//	Cells cells;
	//	auto nb_cells() const {
	//		return cells.nodes.size();
	//	}
	//};

	/** nn stands for number of nodes
	A rotation is performed on the nodes number so that the smallest is the first.
	*/
	template <::std::size_t nn>
	class FaceIndex {
	public:
		typedef ::std::array<NodeId, nn> NodesList;
	private:
		NodesList nodes;
	public:
		FaceIndex() = delete;
		FaceIndex(const NodesList& face_nodes) :
			nodes{ face_nodes } {
			static_assert(nn > 2, "Faces must have 2 nodes at least.");
			auto smallest = ::std::min_element(nodes.begin(), nodes.end());
			::std::rotate(nodes.begin(), smallest, nodes.end());
			// Once the smallest node id comes first there is two way round the face
			// we choose the one that starts with the minimum second id
			if (*::std::next(nodes.begin()) > nodes.back()) {
				::std::reverse(::std::next(nodes.begin()), nodes.end());
			}
		}
		bool operator<(const FaceIndex& other) const {
			return nodes < other.nodes;
		}
		bool operator==(const FaceIndex& other) const {
			return nodes == other.nodes;
		}
		const auto& get_nodes() const {
			return nodes;
		}
	};

	class FaceNeighbors {
	private:
		::std::pair<CellId, CellId> neighbors;
	public:
		FaceNeighbors() = delete;
		FaceNeighbors(const CellId& cell) :
			neighbors{ cell, cell } {}
		bool is_on_boundary() const {
			return neighbors.first == neighbors.second;
		}
		bool is_inside() const {
			return neighbors.first != neighbors.second;
		}
		void associate(const CellId& cell) {
			assert(is_on_boundary());
			neighbors.second = cell;
			assert(is_inside());
		}
		const auto& value() const { return neighbors; }
	};

	template <typename ElementType>
	struct MeshConnectivity {
		typedef ElementType element_type;
		struct Cells {
			FSCoC<NodeId, element_type::nb_nodes()> nodes;
			FSCoC<FaceId, element_type::nb_facets()> faces;
			auto nb() const {
				return nodes.size();
			}
		};
		struct Faces {
			typedef decltype(ElementType::info::facets) Facets;
			typedef typename Facets::value_type Facet;
			static constexpr auto nfn = element_type::nb_facet_nodes();
			typedef FaceIndex<nfn> Index;
			typedef ::std::map<Index, FaceId> Id_factory;
			typedef typename Id_factory::value_type Factory_key;
			Id_factory ids;
			::std::vector<FaceNeighbors> cells;
			FSCoC<NodeId, nfn> nodes;
			template <typename CellNodes>
			void update_from_cellnodes(const CellNodes& cellnodes) {
				ids.clear();
				cells.clear();
				::std::array < NodeId, nfn > facet_nodes;
				const auto ncells = cellnodes.size(); // FIXME: I'd rather use ::std::size(cellnodes) but compilation fails with gcc
				for (CellId cellid = 0; cellid < ncells; ++cellid) {
					for (int facet = 0; facet < element_type::nb_facets(); ++facet) {
						element_type::facet_nodes(cellnodes[cellid], facet, facet_nodes);
						assert(ids.size() == cells.size());
						assert(cells.size()<=std::numeric_limits<FaceId>::max());
						auto face_id = static_cast<FaceId>(cells.size());
						auto key = Factory_key{ Index{ facet_nodes }, face_id };
						auto insertion_result = ids.insert(key);
						auto success = insertion_result.second;
						if (success) {
							cells.emplace_back(cellid);
						}
						else {
							face_id = FaceId{ insertion_result.first->second };
							cells[face_id].associate(cellid);
						}
						assert(ids.size() == cells.size());
					}
				}
				nodes.resize(nb());
				for (auto&& face : ids) {
					nodes[face.second] = face.first.get_nodes();
				}
			}
			auto nb() const {
				return cells.size();
			}
		};
		Cells cells;
		Faces faces;
		void collect_cell_faces() {
			typedef typename Faces::Index FIndex;
			const auto ncells = cells.nb();
			cells.faces.resize(ncells);
			::std::array< NodeId, Faces::nfn > facet_nodes;
			for (CellId cellid = 0; cellid < ncells; ++cellid) {
				for (int facet = 0; facet < element_type::nb_facets(); ++facet) {
					element_type::facet_nodes(cells.nodes[cellid], facet, facet_nodes);
					cells.faces[cellid][facet] = faces.ids.at(FIndex{ facet_nodes });
				}
			}
		}
		void update_from_cellnodes() {
			faces.update_from_cellnodes(cells.nodes);
			collect_cell_faces();
		}
		auto boundary_faces() const {
			::std::vector<FaceId> result;
			const auto nbfaces = faces.nb();
			for (FaceId face = 0; face < nbfaces; ++face) {
				if (faces.cells[face].is_on_boundary()) {
					result.emplace_back(face);
				}
			}
			return result;
		}
		//auto faces_nodes() const {
		//	const auto nbfaces = faces.nb();
		//	::std::vector<::std::array<NodeId, Faces::nfn>> result{ nbfaces };
		//	assert(result.size() == nbfaces);
		//	for (auto&& face : faces.ids) {
		//		result[face.second] = face.first.get_nodes();
		//	}
		//	return result;
		//}
	};

	template <typename VertexType, typename ElementType>
	struct Mesh {
		typedef VertexType Vertex;
		typedef typename Vertex::value_type Coordinate;
		typedef ElementType Element;
		typedef MeshConnectivity<Element> Connectivity;
		::std::vector<Vertex> vertices;
		Connectivity connectivity;
	};

	//typedef LightweightHomogeneousMesh<Tet> TetLMesh;
	//template <sdt::size_t dim>
	//struct Point {
	//	static constexpr auto dimension = dim;
	//	typedef ::std::array<double, dim> Coordinates;
	//	typedef Coordinates::size_type Pos;
	//	Coordinates X;
	//	Point() {
	//		for (Pos i = 0; i < dim; ++i)
	//			X[i] = 0;
	//	}
	//	//Point& operator=(const Point& other) {
	//
	//	//}
	//	static constexpr auto origin() { return Point{}; }
	//
	//
	//};
	typedef ::std::array<double, 3> Point;
	typedef Mesh<Point, Tet> TetMesh;
	typedef Mesh<Point, Hex> HexMesh;

	/** FIXME/IMPROVE Temporary because referenced elements are not copied in memory.
	This might be safer using shared pointers.
	To be compliant with STL algorithm, should define :
	typename _Iter::iterator_category,
	typename _Iter::value_type,
	typename _Iter::difference_type,
	typename _Iter::pointer,
	typename _Iter::reference
	*/
	template <typename Container, typename Indices>
	class TemporarySelection {
	protected:
		Container& container;
		const Indices& indices;
		struct iterator {
			static constexpr auto iterator_category = ::std::forward_iterator_tag{};
			decltype(indices.begin()) pos;
			Container& container;
			iterator() = delete;
			template <typename Position>
			iterator(Container& C, Position&& p) :
				container{ C },
				pos{ ::std::forward<Position>(p) }
			{}
			iterator& operator++() {
				++pos;
				return *this;
			}
			bool operator<(const iterator& other) const {
				assert(&container == &other.container);
				return pos < other.pos;
			}
			bool operator==(const iterator& other) const {
				assert(&container == &other.container);
				return pos == other.pos;
			}
			bool operator!=(const iterator& other) const {
				assert(&container == &other.container);
				return pos != other.pos;
			}
			auto operator*() {
				return container[*pos];
			}
		};
	public:
		TemporarySelection(Container& C, const Indices& I) :
			container{ C },
			indices{ I }
		{}
		auto begin() {
			return iterator{ container, indices.begin() };
		}
		auto end() {
			return iterator{ container, indices.end() };
		}
		auto size() const {
			return indices.size();
		}
	};

	template <typename Container, typename Indices>
	auto make_selection(Container& container, const Indices& indices) {
		return TemporarySelection<Container, Indices>{container, indices};
	}

	template <typename Container, ::std::size_t N, typename IndexType = typename Container::size_type>
	auto extract(Container& container, const ::std::array<IndexType, N>& indices) {
		typedef ::std::array<typename Container::value_type, N> Result;
		Result result;
		for (typename Result::size_type i = 0; i < N; ++i) result[i] = container[indices[i]];
		return result;
	}

	template <typename T, ::std::size_t dim>
	struct Origin {};

	template <typename T>
	struct Origin<T, 2> {
		static constexpr auto build() -> ::std::array<T, 2> {
			return ::std::array<T, 2>{ {0, 0} };
		}
	};

	template <typename T>
	struct Origin<T, 3> {
		static constexpr auto build() -> ::std::array<T, 3> {
			return ::std::array<T, 3>{ {0, 0, 0} };
		}
	};

	template <typename T, ::std::size_t dim>
	constexpr auto origin() -> ::std::array<T, dim>
	{
		return Origin<T, dim>::build();
	}

	template <typename PointIterator>
	auto compute_center(PointIterator first, PointIterator last) {
		constexpr auto dim = ::std::tuple_size<Point>::value;
		static_assert(dim > 0, "Zero dimension for points!");
		double nbpts = 0;
		auto result = Point{
				::std::accumulate(
					first, last, origin<Point::value_type, dim>(),
					[dim, &nbpts](Point& R, const Point& P) {
				for (Point::size_type i = 0; i < dim; ++i) R[i] += P[i];
				++nbpts;
				return R;
		}
				)
		};
		assert(nbpts > 0);
		auto tmp = 1 / nbpts;
		for (Point::size_type i = 0; i < dim; ++i) result[i] *= tmp;
		return result;
	}

	template <typename PointContainer>
	auto compute_center(const PointContainer& points)
	{
		return compute_center(::std::begin(points), ::std::end(points));
	}

	/** Converter from homogeneous table to COC. */
	template <typename T, ::std::size_t N>
	auto FSCoC_as_COC(const FSCoC<T, N>& fscoc)
	{
		std::vector<ElementId> pointers;
		auto n = fscoc.size();
		pointers.reserve(n + 1);
		ElementId pos = 0;
		pointers.emplace_back(pos);
		for (; n != 0; --n) {
			pos += N;
			pointers.emplace_back(pos);
		}
		return std::make_tuple(pointers, fscoc.data()->data());
	}

} // namespace MeshTools

