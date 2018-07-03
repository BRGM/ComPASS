#include <cassert>
#include <vector>
#include <map>
#include <CGAL/Bbox_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

struct Id_info {
    typedef int Id_type;
    static constexpr Id_type default_id = -1;
    Id_type id;
    Id_info() :
        id{ default_id } {}
};

typedef CGAL::Epick Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<Id_info, Kernel> Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Segment Segment;
static_assert(Point::Ambient_dimension::value == 2, "wrong dimension");

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "GenericConnectedComponents.h"

namespace py = pybind11;

template <typename BoxedElement>
struct Segment_intersection_report {
    typedef typename BoxedElement::Handle Element_handle;
    typedef std::multimap < Element_handle, Point > Intersections_map;
    Intersections_map& points;
    Segment_intersection_report(Intersections_map& pts) :
        points{ pts }
    {}
    // callback functor that reports all truly intersecting triangles
    void operator()(const BoxedElement& a, const BoxedElement& b)
    {
        auto intersection = CGAL::intersection(*(a.handle()), *(b.handle()));
        if (intersection) {
            auto P = boost::get<Point>(&(*intersection));
            if (P) {
                points.insert(std::make_pair(a.handle(), *P));
                points.insert(std::make_pair(b.handle(), *P));
            }
            else {
                std::cerr << "CHECK: colinear segments ?????" << std::endl;
            }
        }
    }
};

template <typename BoxedElement, typename Soup>
auto soup_intersection(Soup& soup1, Soup& soup2)
{
    auto put_in_boxes = [](Soup& soup) {
        std::vector<BoxedElement> boxes;
        for (auto p = begin(soup); p != end(soup); ++p)
            boxes.emplace_back(p->bbox(), p);
        return boxes;
    };
    typedef Segment_intersection_report<BoxedElement> Report;
    typedef typename Report::Intersections_map Intersections_map;
    auto points = Intersections_map{};
    auto report = Report{ points };
    auto boxes1 = put_in_boxes(soup1);
    auto boxes2 = put_in_boxes(soup2);
    // Run the self intersection algorithm with all defaults
    CGAL::box_intersection_all_pairs_d(
        begin(boxes1), end(boxes1),
        begin(boxes2), end(boxes2),
        report
    );
    return points;
}

template <typename Input_iterator> 
auto split_segment(Input_iterator p, Input_iterator q)
{
    auto S = p->first;
    auto O = S->source();
    auto u = S->to_vector();
    std::vector<Point> all_points;
    all_points.emplace_back(O);
    std::transform(p, q, std::back_inserter(all_points), [](const typename Input_iterator::value_type& r) { return r.second; });
    all_points.emplace_back(S->target());
    auto less_along_S = [O, u](const Point& P, const Point& Q) { return u*(P - O) < u*(Q - O); };
    std::sort(begin(all_points), end(all_points), less_along_S);
    return all_points;
}

struct Neighborhood
{
    CDT * const cdt;
    CDT::Face_handle face;
    Neighborhood(CDT *p, CDT::Face_handle f) :
        cdt{ p },
        face{ f } {
        assert(p != nullptr);
    }
    static constexpr int size() { return 3; }
    auto operator()(int i) const {
        assert(i < size());
        return face->neighbor(i);
    }
    bool is_valid(int i) const {
        assert(i < size());
        assert(cdt != nullptr);
        return !cdt->is_infinite(face->neighbor(i));
    }
    auto is_connected(int i) const {
        assert(i < size());
        return is_valid(i) && !face->is_constrained(i);
    }
};

struct NeighborhoodFactory
{
    CDT * const cdt;
    typedef CDT::Face_handle Element_handle;
    NeighborhoodFactory(CDT *p) :
        cdt{ p } {
        assert(cdt != nullptr);
    }
    auto operator()(CDT::Face_handle f) const {
        return Neighborhood{ cdt, f };
    }
};

template <typename SegmentSoup>
auto mesh_segment_soup(const SegmentSoup& soup)
{
    CDT cdt;
    for (auto&& S : soup)
        cdt.insert_constraint(S.source(), S.target());

    const std::size_t nv = cdt.number_of_vertices();
    auto vertices = py::array_t<double, py::array::c_style>{ { nv, static_cast<std::size_t>(2) } };
    auto rv = vertices.mutable_unchecked<2>();
    std::map<CDT::Vertex_handle, std::size_t> vmap;
    std::size_t n = 0;
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        vmap[v] = n;
        auto P = v->point();
        for (int k = 0; k<2; ++k) rv(n, k) = P[k];
        ++n;
    }
    const std::size_t nf = cdt.number_of_faces();
    auto triangles = py::array_t<std::size_t, py::array::c_style>{ { nf, static_cast<std::size_t>(3) } };
    auto rt = triangles.mutable_unchecked<2>();
    n = 0;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        for (int k = 0; k<3; ++k) rt(n, k) = vmap[f->vertex(k)];
        ++n;
    }
    typedef ConnectedComponents<NeighborhoodFactory::Element_handle> Connected_components;
    auto components = Connected_components{ 
        cdt.finite_faces_begin(), cdt.finite_faces_end(), NeighborhoodFactory{ &cdt }
    };
    std::map<CDT::Face_handle, std::size_t> fmap;
    n = 0;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        fmap[f] = n;
        ++n;
    }
    auto components_id = py::array_t<std::size_t, py::array::c_style>{ { nf } };
    auto rid = components_id.mutable_unchecked<1>();
    n = 0;
    for (auto&& component : components) {
        for (auto&& f : component) {
            rid(fmap.at(f)) = n;
        }
        ++n;
    }
    return py::make_tuple(vertices, triangles, components_id);
}

auto split_and_mesh(
    py::array_t<double, py::array::c_style> vertices,
    py::array_t<int, py::array::c_style> segments_front,
    py::array_t<int, py::array::c_style> segments_back
) {
    typedef std::vector<Segment> SegmentSoup;
    auto nodes = reinterpret_cast<const Point *>(vertices.unchecked<2>().data(0, 0));
    auto build_segments = [nodes](py::array_t<int, py::array::c_style> segments) {
        auto pairs = segments.unchecked<2>();
        const auto n = pairs.shape(0);
        SegmentSoup v;
        for (int k = 0; k < n; ++k)
            v.emplace_back(*(nodes + pairs(k, 0)), *(nodes + pairs(k, 1)));
        return v;
    };
    auto S1 = build_segments(segments_front);
    auto S2 = build_segments(segments_back);
    typedef CGAL::Box_intersection_d::Box_with_handle_d<
        CGAL::Epick::FT,
        Point::Ambient_dimension::value,
        SegmentSoup::const_iterator> SBox;
    auto intersections = soup_intersection<SBox>(S1, S2);
    auto splitted_segments = SegmentSoup{};
    for (auto p = begin(intersections); p != end(intersections); ++p) {
        auto S = p->first;
        auto q = p; 
        while (q != end(intersections) && q->first == S) ++q;
        auto all_points = split_segment(p, q);
        assert(all_points.size() > 2);
        auto a = begin(all_points);
        for (auto b = next(a); b != end(all_points); ++b) {
            splitted_segments.emplace_back(*a, *b);
            a = b;
        }
    }
    auto add_non_splitted_segments = [&intersections, &splitted_segments](const SegmentSoup& soup) {
        for (auto p = cbegin(soup); p != cend(soup); ++p)
            if (intersections.erase(p) == 0)
                splitted_segments.emplace_back(*p);
    };
    add_non_splitted_segments(S1);
    add_non_splitted_segments(S2);
    return mesh_segment_soup(splitted_segments);
}

auto build_constrained_delaunay_triangulation(
    py::array_t<double, py::array::c_style> vertices,
    py::array_t<int, py::array::c_style> segments
) {
    CDT cdt;
    auto nodes = reinterpret_cast<const Point *>(vertices.unchecked<2>().data(0, 0));
    auto pairs = segments.unchecked<2>();
    const auto n = pairs.shape(0);
    for (int k = 0; k < n; ++k) {
        auto id = pairs(k, 0);
        auto v0 = cdt.insert(*(nodes + id));
        v0->info().id = id;
        id = pairs(k, 1);
        auto v1 = cdt.insert(*(nodes + id));
        v1->info().id = id;
        cdt.insert_constraint(v0, v1);
    }
    return cdt;
}

auto mesh(
    py::array_t<double, py::array::c_style> segment_vertices,
    py::array_t<int, py::array::c_style> segments
) {
    auto cdt = build_constrained_delaunay_triangulation(segment_vertices, segments);
    const std::size_t nv = cdt.number_of_vertices();
    auto vertices = py::array_t<double, py::array::c_style>{ { nv, static_cast<std::size_t>(2) } };
    std::map<CDT::Vertex_handle, std::size_t> vmap;
    std::size_t n = segment_vertices.shape(0);
    auto sp = reinterpret_cast<Point*>(segment_vertices.request().ptr);
    auto p = reinterpret_cast<Point*>(vertices.request().ptr);
    std::copy(sp, sp + n, p);
    assert(n <= nv);
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        auto id = v->info().id;
        if (id < 0) {
            vmap[v] = n;
            *(p + n) = v->point();
            ++n;
        } else {
            vmap[v] = static_cast<std::size_t>(id);
        }
    }
    assert(n == nv);
    const std::size_t nf = cdt.number_of_faces();
    auto triangles = py::array_t<std::size_t, py::array::c_style>{ { nf, static_cast<std::size_t>(3) } };
    auto rt = triangles.mutable_unchecked<2>();
    n = 0;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        for (int k = 0; k<3; ++k) rt(n, k) = vmap[f->vertex(k)];
        ++n;
    }
    typedef ConnectedComponents<NeighborhoodFactory::Element_handle> Connected_components;
    auto components = Connected_components{
        cdt.finite_faces_begin(), cdt.finite_faces_end(), NeighborhoodFactory{ &cdt }
    };
    std::map<CDT::Face_handle, std::size_t> fmap;
    n = 0;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        fmap[f] = n;
        ++n;
    }
    auto components_id = py::array_t<std::size_t, py::array::c_style>{ { nf } };
    auto rid = components_id.mutable_unchecked<1>();
    n = 0;
    for (auto&& component : components) {
        for (auto&& f : component) {
            rid(fmap.at(f)) = n;
        }
        ++n;
    }
    std::set<std::pair<CDT::Vertex_handle, CDT::Vertex_handle>> constrained_edges;
    std::vector<CDT::Vertex_handle> face_nodes;
    py::list all_faces;
    for (auto&& component : components) {
        constrained_edges.clear();
        for (auto&& face : component) {
            for (int i = 0; i < 3; ++i) {
                auto edge = CDT::Edge{ face, i };
                if (cdt.is_constrained(edge) || cdt.is_infinite(face->neighbor(i))) {
                    constrained_edges.insert(
                    { face->vertex((i + 1) % 3), face->vertex((i + 2) % 3) }
                    );
                }
            }
        }
        assert(constrained_edges.size() > 2);
        face_nodes.clear();
        auto edge = begin(constrained_edges);
        face_nodes.push_back(edge->first);
        face_nodes.push_back(edge->second);
        constrained_edges.erase(edge);
        while (!constrained_edges.empty()) {
            for (edge = begin(constrained_edges); edge != end(constrained_edges); ++edge) {
                if (face_nodes.back() == edge->first) {
                    face_nodes.push_back(edge->second);
                    constrained_edges.erase(edge);
                    break;
                }
                if (face_nodes.back() == edge->second) {
                    face_nodes.push_back(edge->first);
                    constrained_edges.erase(edge);
                    break;
                }
            }
        }
        assert(face_nodes.front() == face_nodes.back());
        face_nodes.pop_back();
        py::list l;
        for (auto&& v : face_nodes) {
            l.append(vmap[v]);
        }
        all_faces.append(l);
    }
    return py::make_tuple(vertices, triangles, components_id, all_faces);
}

PYBIND11_MODULE(PetrelMesh, module)
{

    module.doc() = "pybind11 with CGAL2D (quick and dirty!!!)";
    module.def("split_and_mesh", &split_and_mesh);
    module.def("mesh", &mesh);
    //module.def("foo", []() {
    //    boost::optional<int> test;
    //    test = false;
    //    if (test) py::print("ya bon");
    //});
}
