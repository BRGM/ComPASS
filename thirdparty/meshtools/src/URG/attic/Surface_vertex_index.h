#pragma once

template <typename Surface_index>
struct Surface_vertex_index
{
    // There would be much less lines of code direcly replacing the variant by a std::set<Surface_index>
    //struct Not_indexed {};
    //typedef boost::variant <
    //    Not_indexed,
    //    Surface_index,
    //    std::array<Surface_index, 2>,
    //    std::array<Surface_index, 3>,
    //    std::set<Surface_index>
    //> Index_type;
    //Index_type index;
    //struct Adder : boost::static_visitor<Index_type> {
    //    Surface_index added;
    //    Adder(Surface_index si) :
    //        boost::static_visitor<Index_type>{},
    //        added{ si }
    //    {}
    //    result_type operator()(const Not_indexed&) const {
    //        return added;
    //    }
    //    result_type operator()(const Surface_index& si) const {
    //        if (added == si) return si;
    //        if (si < added) return std::array<Surface_index, 2>{ {si, added } };
    //        return std::array<Surface_index, 2>{ {added, si } };
    //    }
    //    result_type operator()(const std::array<Surface_index, 2>& sia) const {
    //        assert(sia[0] < sia[1]);
    //        if (added == sia[0] || added == sia[1]) return sia;
    //        if (added < sia[0]) return std::array<Surface_index, 3>{ { added, sia[0], sia[1] } };
    //        if (added < sia[1]) return std::array<Surface_index, 3>{ { sia[0], added, sia[1] } };
    //        return std::array<Surface_index, 3>{ { sia[0], sia[1], added } };
    //    }
    //    result_type operator()(const std::array<Surface_index, 3>& sia) const {
    //        if (added == sia[0] || added == sia[1] || added == sia[2]) return sia;
    //        auto s = std::set<Surface_index>{ begin(sia), end(sia) };
    //        s.insert(added);
    //        return s;
    //    }
    //    result_type operator()(std::set<Surface_index>& s) const {
    //        s.insert(added);
    //        return s;
    //    }
    //};
    //struct Degree : boost::static_visitor<std::size_t> {
    //    constexpr result_type operator()(const Not_indexed&) const noexcept {
    //        return 0;
    //    }
    //    constexpr result_type operator()(const Surface_index&) const noexcept {
    //        return 1;
    //    }
    //    constexpr result_type operator()(const std::array<Surface_index, 2>&) const noexcept {
    //        return 2;
    //    }
    //    constexpr result_type operator()(const std::array<Surface_index, 3>&) const noexcept {
    //        return 3;
    //    }
    //    result_type operator()(std::set<Surface_index>& s) const noexcept {
    //        return s.size();
    //    }
    //};
    //Surface_vertex_index& add(Surface_index si) {
    //    index = boost::apply_visitor(Adder{ si }, index);
    //    return *this;
    //}
    //std::size_t degree() const {
    //    return boost::apply_visitor(Degree{}, index);
    //}
    typedef std::set<Surface_index> Index_type;
    Index_type index;
    Surface_vertex_index() = default;
    Surface_vertex_index(const Surface_vertex_index&) = default;
    Surface_vertex_index& operator=(const Surface_vertex_index&) = default;
    Surface_vertex_index(Surface_index si) :
        index{ si }
    {}
    Surface_vertex_index& add(Surface_index si) {
        index.insert(si);
        return *this;
    }
    Surface_vertex_index& remove(Surface_index si) {
        index.erase(si);
        return *this;
    }
    auto degree() const {
        return index.size();
    }
};

template <typename Surface_index>
bool operator==(const Surface_vertex_index<Surface_index>& svi1, const Surface_vertex_index<Surface_index>& svi2)
{
    const auto& s1 = svi1.index;
    const auto& s2 = svi2.index;
    if (s1.size() != s2.size()) return false;
    std::vector<Surface_index> intersection;
    intersection.reserve(s1.size());
    std::set_intersection(begin(s1), end(s1), begin(s2), end(s2), std::back_inserter(intersection));
    return intersection.size() == s1.size();
}

//template <typename Surface_index, typename Surface_mesh>
//auto add_vertex_degree_map(Surface_mesh& mesh)
//{
//    typedef typename Surface_mesh::Vertex_index Vertex_index;
//    typedef Surface_vertex_index<Surface_index> Vertex_degree;
//    auto pmap = mesh.add_property_map<Vertex_index, Vertex_degree>("vertex_degree", Vertex_degree{});
//    assert(pmap.second);
//    return pmap.first;
//}
//
//template <typename Surface_mesh>
//auto add_is_contrained_edge_map(Surface_mesh& mesh)
//{
//    typedef typename Surface_mesh::Edge_index Edge_index;
//    auto pmap = mesh.add_property_map<Edge_index, bool>("is_constrained_edge", false);
//    assert(pmap.second);
//    return pmap.first;
//}

//template <typename Surface_index, typename Surface_mesh>
//void add_properties(Surface_mesh& mesh)
//{
//    add_vertex_degree_map<Surface_index>(mesh);
//    add_is_contrained_edge_map(mesh);
//}

