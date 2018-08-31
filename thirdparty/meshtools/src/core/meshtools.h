#pragma once

#include <cassert>
#include <array>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <iterator>
#include <set>

#include "common.h"
#include "Variant_utils.h"

/* TODO ?
List facets as vector instead of map (iteration over faces).
Might recopy in algorithm ?
We could sort faces according to Face_index (same order as in map) but
this would make addition / mesh merging more difficult (order might change)
*/

namespace MeshTools
{

	typedef ElementId NodeId;
	typedef ElementId CellId;
	typedef ElementId FaceId;
	// FIXME: We'd better use std::size_t here to avoid static casts
	//        this is for compatibility reasons with ComPASS internal Fortran mesh structures
	typedef ElementId DefaultCOCPointerType;
	
	/** Elements Connectivity is exactly the same as in the VTK library.
	cf. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf */
	// FIXME: Use enum instead of VTK id???
	typedef int VTK_ID_type;

	template <std::size_t nn>
	auto change_nodes(const std::array<NodeId, nn>& input,
					  const std::map<NodeId, NodeId>& new_ids)
	{
		auto output = std::array<NodeId, nn>{};
		for (std::size_t k = 0; k < nn; ++k) {
			output[k] = new_ids.at(input[k]);
		};
		return output;
	}

	template <std::size_t nn>
	struct Element_by_nodes
	{
		std::array<NodeId, nn> nodes;
		// CHECKME: We fix the first constructor parameter not to shadow other constructor (copy constructor...)
		//          Otherwise variant compilation fails with gcc
		template <typename ... NodeIds>
		Element_by_nodes(NodeId n0, NodeIds... nis) :
			nodes{ { n0, nis... } }
		{
			static_assert(sizeof...(NodeIds) == nn - 1, "Wrong nomber of nodes!");
			static_assert(is_single_type<std::tuple<NodeIds...>>(), "Wrong type for nodes ids!");
			static_assert(std::is_same<NodeId, typename IsSingleTypeTuple<std::tuple<NodeIds...>>::First_element_type>::value, "Wrong type for nodes ids!");
		}
		Element_by_nodes() :
			nodes{}
		{}
		Element_by_nodes(const std::array<NodeId, nn>& nodes_array) :
			nodes{ nodes_array }
		{}
		Element_by_nodes(std::array<NodeId, nn>&& nodes_array) :
			nodes{ std::forward<std::array<NodeId, nn>>(nodes_array) }
		{}
		constexpr static std::size_t nbnodes() { return nn; }
		constexpr bool operator==(const Element_by_nodes& other) const noexcept {
			return nodes == other.nodes;
		}
	};

	template <typename FacetTypes>
	using FacetsType = VariantFromTupleIfNeeded<FacetTypes>;

	template <typename ... Facets>
	struct With_facets
	{
		typedef std::tuple<Facets...> Facet_types;
		typedef FacetsType<Facet_types> Facets_type;
		typedef std::array<Facets_type, std::tuple_size<Facet_types>::value> Facets_array;
		constexpr static std::size_t nbfacets() { return std::tuple_size<Facet_types>::value; }
	};

    struct Line :
        Element_by_nodes<2>
    {
        using Element_by_nodes<nbnodes()>::Element_by_nodes;
        static constexpr auto name = "Line";
        static constexpr auto VTK_ID = VTK_ID_type{ 3 };
    };

    struct Triangle :
        Element_by_nodes<3>,
        With_facets<Line, Line, Line>
    {
        using Element_by_nodes<nbnodes()>::Element_by_nodes;
        static constexpr auto name = "Triangle";
        static constexpr auto VTK_ID = VTK_ID_type{ 5 };
        auto facets() const -> Facets_array {
            return Facets_array{
                Line{ nodes[0], nodes[1] },
                Line{ nodes[1], nodes[2] },
                Line{ nodes[2], nodes[0] }
            };
        };
    };

    struct Quad :
		Element_by_nodes<4>
	{
		using Element_by_nodes<nbnodes()>::Element_by_nodes;
		static constexpr auto name = "Quad";
		static constexpr auto VTK_ID = VTK_ID_type{ 9 };
	};

	struct Tetrahedron :
		Element_by_nodes<4>,
		With_facets<Triangle, Triangle, Triangle, Triangle>
	{
		using Element_by_nodes<nbnodes()>::Element_by_nodes;
		static constexpr auto name = "Tetrahedron";
		static constexpr auto VTK_ID = VTK_ID_type{ 10 };
		// Outside oriented facets
		auto facets() const -> Facets_array {
			return Facets_array{
				Triangle{ nodes[1], nodes[2], nodes[3] },
				Triangle{ nodes[0], nodes[3], nodes[2] },
				Triangle{ nodes[0], nodes[1], nodes[3] },
				Triangle{ nodes[0], nodes[2], nodes[1] }
			};
		};
	};

	struct Hexahedron :
		Element_by_nodes<8>,
		With_facets<Quad, Quad, Quad, Quad, Quad, Quad>
	{
		using Element_by_nodes<nbnodes()>::Element_by_nodes;
		static constexpr auto name = "Hexahedron";
		static constexpr auto VTK_ID = VTK_ID_type{ 12 };
		// Outside oriented facets
		auto facets() const -> Facets_array {
			return Facets_array{
				Quad{ nodes[0], nodes[1], nodes[2], nodes[3] },
				Quad{ nodes[4], nodes[5], nodes[6], nodes[7] },
				Quad{ nodes[1], nodes[2], nodes[6], nodes[5] },
				Quad{ nodes[2], nodes[6], nodes[7], nodes[3] },
				Quad{ nodes[3], nodes[7], nodes[4], nodes[0] },
				Quad{ nodes[0], nodes[1], nodes[5], nodes[4] }
			};
		};
	};

    struct Wedge :
        Element_by_nodes<6>,
        With_facets<Triangle, Triangle, Quad, Quad, Quad>
    {
        using Element_by_nodes<nbnodes()>::Element_by_nodes;
        static constexpr auto name = "Wedge";
        static constexpr auto VTK_ID = VTK_ID_type{ 13 };
        // Outside oriented facets
        auto facets() const -> Facets_array {
            return Facets_array{
                Triangle{ nodes[0], nodes[1], nodes[2] },
                Triangle{ nodes[3], nodes[4], nodes[5] },
                Quad{ nodes[0], nodes[2], nodes[5], nodes[3] },
                Quad{ nodes[1], nodes[4], nodes[5], nodes[2] },
                Quad{ nodes[0], nodes[3], nodes[4], nodes[1] }
            };
        };
    };

    struct Pyramid :
        Element_by_nodes<5>,
        With_facets<Quad, Triangle, Triangle, Triangle, Triangle>
    {
        using Element_by_nodes<nbnodes()>::Element_by_nodes;
        static constexpr auto name = "Pyramid";
        static constexpr auto VTK_ID = VTK_ID_type{ 14 };
        // Outside oriented facets
        auto facets() const -> Facets_array {
            return Facets_array{
                Quad{ nodes[0], nodes[1], nodes[2], nodes[3] },
                Triangle{ nodes[0], nodes[1], nodes[4] },
                Triangle{ nodes[1], nodes[2], nodes[4] },
                Triangle{ nodes[2], nodes[3], nodes[4] },
                Triangle{ nodes[3], nodes[0], nodes[4] },
            };
        };
    };

    template <typename Element>
	constexpr auto change_element_nodes(const Element& element,
		                                   const std::map<NodeId, NodeId>& new_nodes)
	{
		return Element{ change_nodes(element.nodes, new_nodes) };
	}

	template <typename ... Alternatives>
	constexpr auto change_element_nodes(const mapbox::util::variant<Alternatives...>& variant,
		                                const std::map<NodeId, NodeId>& new_nodes)
	{
		typedef mapbox::util::variant<Alternatives...> Element;
		return mapbox::util::apply_visitor(
			[&new_nodes](auto alternative) {
				return Element{
					change_element_nodes(alternative, new_nodes)
				};
		    },
		    variant);
	}

	template <typename FacetType, typename ElementType>
	inline auto facet_index(const FacetType& facet, const ElementType& element)
	{
		auto facets = element.facets();
		for (std::size_t k = 0; k!=ElementType::nbfacets(); ++k) {
			if (facet == facets[k]) return k;
		}
		assert(false);
		return ElementType::nbfacets();
	}

	struct ArrayBaseAccessor
	{
		template <typename T, std::size_t n>
		static constexpr std::size_t size(const std::array<T, n>&) noexcept
		{
			return n;
		}
		template <typename T, std::size_t n>
		static constexpr const T * begin(const std::array<T, n>& a) noexcept
		{
			return a.data();
		}
		template <typename T, std::size_t n>
		static constexpr T& access(std::array<T, n>& a, std::size_t k) noexcept
		{
			assert(k < n);
			return a[k];
		}
		template <typename T, std::size_t n>
		static constexpr T extract(const std::array<T, n>& a, std::size_t k) noexcept
		{
			assert(k < n);
			return a[k];
		}
	};

	// FIXME/IMPROVE: This is necessary only because Element_by_nodes does not inherit from std::array<NodeId, n> 
	struct NodesBaseAccessor
	{
		template <typename Element>
		static constexpr std::size_t size(const Element&) noexcept
		{
			return Element::nbnodes();
		}
		template <typename Element>
		static constexpr const NodeId * begin(const Element& element) noexcept
		{
			return element.nodes.data();
		}
		template <typename Element>
		static constexpr NodeId& access(Element& element, std::size_t k) noexcept
		{
			assert(k < Element::nbnodes());
			return element.nodes[k];
		}
		template <typename Element>
		static constexpr NodeId extract(const Element& element, std::size_t k) noexcept
		{
			assert(k < Element::nbnodes());
			return element.nodes[k];
		}
	};

	template <typename BaseAccessor>
	struct VariantAccessor : BaseAccessor
	{
		using BaseAccessor::size;
		template <typename ... Alternatives>
		static constexpr auto size(const mapbox::util::variant<Alternatives...>& variant) noexcept
		{
			return mapbox::util::apply_visitor([](auto alternative) { return BaseAccessor::size(alternative); }, variant);
		}
		using BaseAccessor::begin;
		template <typename ... Alternatives>
		static constexpr auto begin(const mapbox::util::variant<Alternatives...>& variant) noexcept
		{
			// we need reference here because we return a pointer and we don't want to make a copy of element
			return mapbox::util::apply_visitor([](const auto& alternative) { return BaseAccessor::begin(alternative); }, variant);
		}
		using BaseAccessor::access;
		template <typename ... Alternatives>
		static constexpr auto access(mapbox::util::variant<Alternatives...>& variant, std::size_t k) noexcept
		{
			return mapbox::util::apply_visitor([k](auto& alternative) { return BaseAccessor::access(alternative, k); }, variant);
		}
		using BaseAccessor::extract;
		template <typename ... Alternatives>
		static constexpr auto extract(const mapbox::util::variant<Alternatives...>& variant, std::size_t k) noexcept
		{
			return mapbox::util::apply_visitor([k](const auto& alternative) { return BaseAccessor::extract(alternative, k); }, variant);
		}
	};

	template <typename access_type>
	struct Range
	{
	protected:
		std::pair<access_type, access_type> bounds;
	public:
		Range(access_type first, access_type last) :
			bounds{ first, last }
		{}
		Range(access_type&& first, access_type&& last) :
			bounds{ std::forward<access_type>(first), std::forward<access_type>(last) }
		{}
		auto begin() const { return bounds.first; }
		auto end() const { return bounds.second; }
		auto size() const { return bounds.second - bounds.first; }
	};

	template <typename BaseAccessor = ArrayBaseAccessor, typename Element>
	inline auto range(const Element& element) {
		typedef VariantAccessor<BaseAccessor> Accessor;
		auto const p = Accessor::begin(element);
		return Range<decltype(p)>{ p, p + Accessor::size(element) };
	}

	// FIXME/IMPROVE: This is necessary only because Element_by_nodes does not inherit from std::array<NodeId, n> 
	template <typename Element>
	inline auto node_range(const Element& element) {
		return range<NodesBaseAccessor>(element);
	}

	template <typename Element>
	inline constexpr VTK_ID_type vtk_id(const Element&) noexcept
	{
		return Element::VTK_ID;
	}
	
	template <typename ... Alternatives>
	inline constexpr VTK_ID_type vtk_id(const mapbox::util::variant<Alternatives...>& variant) noexcept
	{
		return mapbox::util::apply_visitor([](const auto& alternative) { return vtk_id(alternative); }, variant);
	}

	template <typename Accessor, typename IDType, typename PointerType = DefaultCOCPointerType, typename Collection>
	inline auto collection_as_COC_data(const Collection& collection) ->
		std::tuple<std::vector<PointerType>, std::vector<IDType>>
	{
		const auto n = collection.size();
		auto pointers = std::vector<PointerType>{};
		pointers.reserve(n + 1);
		PointerType pos = 0;
		pointers.emplace_back(pos);
		for (auto&& element : collection) {
			pos += static_cast<PointerType>(Accessor::size(element));
			pointers.emplace_back(pos);
		}
		assert(pointers.size() == n + 1);
		auto nodes = std::vector<IDType>{};
		nodes.reserve(pos);
		for (auto&& element : collection) {
			auto p = Accessor::begin(element);
			for (std::size_t k = 0; k < Accessor::size(element); ++k) {
				nodes.emplace_back(*p);
				++p;
			}
		}
		return std::make_tuple(pointers, nodes);
	}

	template <typename IDType = NodeId, typename PointerType = DefaultCOCPointerType, typename Collection>
	inline auto node_collection_as_COC_data(const Collection& collection)
	{
		return collection_as_COC_data<VariantAccessor<NodesBaseAccessor>, IDType, PointerType>(collection);
	}

	template <typename RawFacetType>
	struct RawFaceIndex 
	{
		typedef RawFacetType Facet_type;
		typedef decltype(Facet_type::nodes) Index;
	protected:
		Index index;
		/** A rotation is performed on the nodes number so that the smallest is the first. */
		void sort_nodes_for_index()
		{
			static_assert(std::tuple_size<Index>::value, "Faces must have 2 nodes at least.");
			auto smallest = ::std::min_element(index.begin(), index.end());
			::std::rotate(index.begin(), smallest, index.end());
			// Once the smallest node id comes first there is two way round the face
			// we choose the one that starts with the minimum second id
			if (*::std::next(index.begin()) > index.back()) {
				::std::reverse(::std::next(index.begin()), index.end());
			}
		}
	public:
		RawFaceIndex() = delete;
		RawFaceIndex(const RawFacetType& facet) :
			index{ facet.nodes } 
		{
			sort_nodes_for_index();
		}
		auto as_facet() const {
			return Facet_type{ index };
		}
		bool operator<(const RawFaceIndex& other) const {
			return index < other.index;
		}
		bool operator==(const RawFaceIndex& other) const {
			return index == other.index;
		}
	};

	template <typename ... CellTypes>
	using AllRawFacetTypes = typename ConcatenatedTupleSet<typename CellTypes::Facet_types...>::type;

	template <typename RawFacetTuple>
	struct RawFaceIndexTuple;

	template <typename ... RawFacetTypes>
	struct RawFaceIndexTuple < std::tuple<RawFacetTypes...>>
	{	
		typedef std::tuple< RawFaceIndex<RawFacetTypes>... > type;
	};

	template <typename ... CellTypes>
	using FaceIndex = VariantFromTupleIfNeeded<typename RawFaceIndexTuple<AllRawFacetTypes<CellTypes...>>::type>;

	class FaceNeighbors {
	protected:
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

	template <typename CellType>
	using CellFacesId = std::array<FaceId, CellType::nbfacets()>;

	template <typename IDType = FaceId, typename PointerType = DefaultCOCPointerType, typename Collection>
	inline auto face_collection_as_COC_data(const Collection& collection)
	{
		return collection_as_COC_data<VariantAccessor<ArrayBaseAccessor>, IDType, PointerType>(collection);
	}

	template <typename ... CellTypes>
	struct AllCellFacesId;

	template <typename CellType, typename ... CellTypes>
	struct AllCellFacesId<CellType,CellTypes...>
	{
		typedef typename AllCellFacesId<CellTypes...>::type Cell_faces_id_subtuple;
		typedef CellFacesId<CellType> First_element_type;
		typedef typename ExtendTupleSet<Cell_faces_id_subtuple, First_element_type, !IsInTuple<First_element_type, Cell_faces_id_subtuple>::value>::Extended_tuple_type type;
	};

	template <typename CellType>
	struct AllCellFacesId<CellType>
	{
		typedef std::tuple<CellFacesId<CellType>> type;
	};

	template <typename ... CellTypes>
	inline constexpr auto all_facets_number() { 
		return ::std::make_tuple(CellTypes::nbfacets()...); 
	}

	template <typename NewFacetType, typename ... CellTypes>
	struct AllNewFacetArray;

	template <typename NewFacetType, typename CellType, typename ... CellTypes>
	struct AllNewFacetArray<NewFacetType, CellType, CellTypes...>
	{
		typedef typename AllNewFacetArray<NewFacetType, CellTypes...>::type New_facet_array_subtuple;
		typedef ::std::array<NewFacetType, CellType::nbfacets()> First_element_type;
		typedef typename ExtendTupleSet<New_facet_array_subtuple, First_element_type, !IsInTuple<First_element_type, New_facet_array_subtuple>::value>::Extended_tuple_type type;
	};

	template <typename NewFacetType, typename CellType>
	struct AllNewFacetArray<NewFacetType, CellType>
	{
		typedef ::std::tuple<::std::array<NewFacetType, CellType::nbfacets()>> type;
	};

	//template <typename NewFacetType, typename ... CellTypes>
	//using MeshFacetArray = VariantFromTupleIfNeeded<typename AllNewFacetArray<NewFacetType, CellTypes...>::type>;

	template <typename ... CellTypes>
	struct MeshTraits
	{
		typedef std::tuple<CellTypes...> Cell_types;
		typedef VariantIfNeeded<CellTypes...> Cell_type;
		//FIXME: Why the following does not work???
		//typedef VariantIfNeeded<CellFacesId<CellTypes>...> Cell_faces_id;
		typedef VariantFromTupleIfNeeded<typename AllCellFacesId<CellTypes...>::type> Cell_faces_id;
		typedef VariantFromTupleIfNeeded<AllRawFacetTypes<CellTypes...>> Facet_type;
		typedef VariantFromTupleIfNeeded<typename AllNewFacetArray<Facet_type, CellTypes...>::type> Facets_array;
		typedef FaceIndex<CellTypes...> Face_index;
	};
	
	template <typename ... CellTypes>
	struct MeshConnectivity 
	{
		typedef MeshTraits<CellTypes...> Mesh_traits;
		typedef typename Mesh_traits::Cell_type Cell_type;
		template <typename FacetType>
		static auto recast_facet(const FacetType& facet)
		{
			typedef typename Mesh_traits::Facet_type Facet_type;
			return Facet_type{ facet };
		}
		template <typename ... FacetTypes>
		static auto recast_facet(const mapbox::util::variant<FacetTypes...>& vfacet)
		{
			typedef typename Mesh_traits::Facet_type Facet_type;
			return mapbox::util::apply_visitor([](const auto& facet) { return Facet_type{ facet }; }, vfacet);
		}
		template <typename FacetType>
		static auto make_face_index(const FacetType& facet)
		{
			typedef typename Mesh_traits::Face_index Face_index;
			return Face_index{ facet };
		}
		template <typename ... FacetTypes>
		static auto make_face_index(const mapbox::util::variant<FacetTypes...>& vfacet)
		{
			typedef typename Mesh_traits::Face_index Face_index;
			return mapbox::util::apply_visitor([](const auto& facet) { return Face_index{ facet }; }, vfacet);
		}
		// OPTIMIZE: Would it be possible to return a fixed size array
		//           using specialization and a function to register facets
		template <typename CellType>
		static auto extract_cell_facets(const CellType& cell)
		{
			typedef typename Mesh_traits::Facet_type Facet_type;
			auto facets = cell.facets();
			constexpr ::std::size_t n = CellType::nbfacets();
			auto result = ::std::vector<Facet_type>{};
			result.reserve(n);
			// Loop shall be unrolled as n is small and known at compile time (constexpr)
			for (::std::size_t i = 0; i < n; ++i) {
				result.emplace_back(recast_facet(facets[i]));
			}
			return result;
		}
		template <typename ... CellTs>
		static auto extract_cell_facets(const mapbox::util::variant<CellTs...>& vcell)
		{
			return mapbox::util::apply_visitor([](const auto& cell) {
				return extract_cell_facets(cell);
			}, vcell);
		}
		struct Cells 
		{
			typedef typename Mesh_traits::Cell_type Cell_type;
			typedef typename Mesh_traits::Cell_faces_id Cell_faces_id;
			typedef ::std::vector<Cell_type> Nodes;
			typedef ::std::vector<Cell_faces_id> Faces;
			Nodes nodes;
			Faces faces;
			auto nb() const {
				return nodes.size();
			}
		};
		struct Faces
		{
			// We use a map to provide binary search on Face_index
			typedef typename Mesh_traits::Facet_type Facet_type;
			typedef typename Mesh_traits::Face_index Face_index;
			typedef ::std::map<Face_index, FaceId> Id_factory;
			typedef typename Id_factory::value_type Factory_key;
			typedef ::std::vector<FaceNeighbors> Cells;
			typedef ::std::vector<Facet_type> Nodes;
			Id_factory ids;
			Cells cells;
			Nodes nodes;
			template <typename FaceIndexType>
			static auto extract_facet_from_index(const FaceIndexType& index) {
				return Facet_type{ index.as_facet() };
			}
			template <typename ... FaceIndexTypes>
			static auto extract_facet_from_index(const mapbox::util::variant<FaceIndexTypes...>& vindex) {
				return mapbox::util::apply_visitor([](const auto& index) { return extract_facet_from_index(index); }, vindex);
			}
			//template <std::size_t n>
			//void register_facets(const CellId& cellid, const std::array<Facet_type, n>& facets) {
			//}
			//template <typename ... FacetsArray>
			//void register_facets(const CellId& cellid, const mapbox::util::variant<FacetsArray...>& vfacets) {
			//	mapbox::util::apply_visitor([this, &cellid](const auto& facets) { register_facets(cellid, facets); }, vfacets);
			//}
			template <typename Cells>
			void update_from_cellnodes(const Cells& new_cellnodes) {
				ids.clear();
				cells.clear();
				auto cellid = CellId{ 0 };
				for (auto&& cell : new_cellnodes) {
					for (auto&& facet : extract_cell_facets(cell)) {
						assert(ids.size() == cells.size());
						assert(cells.size() <= std::numeric_limits<FaceId>::max());
						auto face_id = static_cast<FaceId>(cells.size());
						auto key = Factory_key{ make_face_index(facet), face_id };
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
					++cellid;
				}
				// We need to sort face indexes according to face id
				::std::map<FaceId, Face_index> indexes;
				for (auto&& face: ids) {
					indexes.emplace(::std::make_pair(face.second, face.first));
				}
				nodes.clear();
				nodes.reserve(nb());
				assert(indexes.size() == nb());
				for (auto&& face: indexes) {
					assert(face.first == nodes.size());
					nodes.emplace_back(extract_facet_from_index(face.second));
				}
			}
			auto nb() const {
				assert(ids.size() == cells.size());
				return cells.size();
			}
		};
		Cells cells;
		Faces faces;
		template <typename CellType>
		auto collect_cell_faces_id(const CellType& cell) const {
			typedef typename Mesh_traits::Cell_faces_id Cell_faces_id;
			std::array<FaceId, CellType::nbfacets()> result;
			std::size_t k = 0;
			for (auto&& facet : cell.facets()) {
				result[k] = faces.ids.at(make_face_index(facet));
				++k;
			}
			assert(k == CellType::nbfacets());
			return Cell_faces_id{ result };
		}
		template <typename ... CellTs>
		auto collect_cell_faces_id(mapbox::util::variant<CellTs...>& vcell) const {
			return mapbox::util::apply_visitor([this](const auto& cell) { return this->collect_cell_faces_id(cell); }, vcell);
		}
		void collect_cell_faces() {
			cells.faces.clear();
			cells.faces.reserve(cells.nb());
			for (auto&& cell: cells.nodes) {
				cells.faces.emplace_back(collect_cell_faces_id(cell));
			}
		}
		auto boundary_cells() const {
			::std::vector<CellId> result;
			auto fk = FaceId{ 0 };
			for (auto&& face : faces.cells) {
				if (face.is_on_boundary()) {
					result.emplace_back(face.value().first);
				}
				++fk;
			}
			assert(fk == faces.nb());
			return result;
		}
		auto boundary_faces() const {
			::std::vector<FaceId> result;
			auto fk = FaceId{ 0 };
			for (auto&& face : faces.cells) {
				if (face.is_on_boundary()) {
					result.emplace_back(fk);
				}
				++fk;
			}
			assert(fk == faces.nb());
			return result;
		}
		void update_from_cellnodes() {
			faces.update_from_cellnodes(cells.nodes);
			collect_cell_faces();
		}
		void remap_cellnodes(
			const Cell_type* first, const Cell_type* last,
			std::map<NodeId, NodeId>& new_nodes
		) {
			cells.nodes.clear();
			cells.faces.clear();
			cells.nodes.reserve(std::distance(first, last));
			std::transform(
				first, last,
				std::back_inserter(cells.nodes),
				[&new_nodes](const Cell_type& cell) {
					return change_element_nodes(cell, new_nodes);
				}
			);
			update_from_cellnodes();
		}
	};

	template <typename Coordinate, std::size_t dim>
	struct Point: ::std::array<Coordinate, dim>
	{
		typedef std::array<Coordinate, dim> Base_array;
		typedef Coordinate Coordinate_type;
		constexpr static auto dimension = dim;
		template <typename ... Coordinates>
		// CHECKME: We fix the first constructor parameter not to shadow other constructor (copy constructor...)
		//          Otherwise variant compilation fails with gcc
		Point(Coordinate_type x0, Coordinates... xis) :
			::std::array<Coordinate, dim>{ { x0, xis... } }
		{
			static_assert(sizeof...(Coordinates) == dim -1, "Wrong nomber of coordinates!");
			static_assert(is_single_type<std::tuple<Coordinates...>>(), "Wrong type for nodes ids!");
			static_assert(std::is_same<Coordinate_type, typename IsSingleTypeTuple<std::tuple<Coordinates...>>::First_element_type>::value, "Wrong type for nodes ids!");
		}
		Point() :
			Base_array{}
		{}
		Point(const Base_array& a) :
			Base_array{ a }
		{}
		Point(Base_array&& a) :
			Base_array{ std::forward<Base_array>(a) }
		{}
		Point& operator=(const Point& other) {
			Base_array::operator=(other);
			return *this;
		}
		Point& operator+=(const Point& other) {
			for (std::size_t i = 0; i < dim; ++i) {
				(*this)[i] += other[i];
			}
			return *this;
		}
		static Point origin() {
			Point result;
			for (std::size_t i = 0; i < dim; ++i) {
				result[i] = 0;
			}
			return result;
		};
	};

	template <typename Coordinate, std::size_t dim>
	struct PointAccumulator
	{
		typedef Point<Coordinate, dim> Point_type;
		Point_type accumulation;
		std::size_t n;
		PointAccumulator() :
			accumulation{ Point_type::origin() },
			n{ 0 }
		{}
		PointAccumulator(const Point_type& P) :
			accumulation{ P },
			n{ 1 }
		{}
		PointAccumulator(Point_type&& P) :
			accumulation{ std::forward<Point_type>(P) },
			n{ 1 }
		{}
		void operator()(const Point_type& P) {
			accumulation += P;
			++n;
		}
		auto average() const {
			assert(n > 0);
			auto div = 1. / static_cast<double>(n);
			auto result = accumulation;
			for (std::size_t i = 0; i < dim; ++i) {
				result[i]*= div;
			}
			return result;
		}
	};

	template <typename VertexType, typename ... CellTypes>
	struct Mesh {
		typedef VertexType Vertex_type;
		typedef typename Vertex_type::value_type Coordinate_type;
		typedef ::std::vector<Vertex_type> Vertices;
		constexpr static auto dimension = Vertex_type::dimension;
		typedef MeshConnectivity<CellTypes...> Connectivity;
		typedef typename Connectivity::Mesh_traits Mesh_traits;
        typedef typename Mesh_traits::Cell_type Cell_type;
        typedef typename Mesh_traits::Facet_type Facet_type;
        Vertices vertices;
		Connectivity connectivity;
		template <typename Functor, typename NodeElementType>
		void apply_on_vertices(Functor& F, const NodeElementType& element) const {
			for (auto&& node: element.nodes) {
				F(vertices[node]);
			}
		}
		template <typename Functor, typename ... NodeElementTypes>
		void apply_on_vertices(Functor& F, const mapbox::util::variant<NodeElementTypes...>& velement) const {
			mapbox::util::apply_visitor([this, &F](const auto& element) { this->apply_on_vertices(F, element); }, velement);
		}
		template <typename ElementContainer>
		auto compute_centers(const ElementContainer& elements) const {
			auto result = Vertices{};
			result.reserve(elements.size());
			for (auto&& element : elements) {
				PointAccumulator<Coordinate_type, dimension> S;
				apply_on_vertices(S, element);
				result.emplace_back(S.average());
			}
			return result;
		}
		auto cell_centers() const {
			return compute_centers(connectivity.cells.nodes);
		}
		auto face_centers() const {
			return compute_centers(connectivity.faces.nodes);
		}
		template <typename MeshElement>
		auto element_vertices(const MeshElement& element) const {
			static_assert(IsInTuple<MeshElement, typename Mesh_traits::Cell_types>::value, "Mesh does not hold elements of this kind!");
			auto result = std::vector<Vertex_type>{};
			result.reserve(element.nbnodes());
			for (auto&& ei : element.nodes) {
				result.emplace_back(vertices[ei]);
			}
			return result;
		}
		auto nb_vertices() const { return vertices.size(); }
		auto nb_cells() const { return connectivity.cells.nb(); }
		auto nb_faces() const { return connectivity.faces.nb(); }
	};

	typedef Mesh<Point<double, 3>, Tetrahedron> TetMesh;
	typedef Mesh<Point<double, 3>, Hexahedron> HexMesh;

} // namespace MeshTools

