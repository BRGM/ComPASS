#pragma once

template <typename Constraint_map>
auto collect_consrained_edges_as_curves(const Mesh& mesh, const Constraint_map& constraints)
{

    typedef typename Mesh::Vertex_index Vertex_index;
    typedef typename Mesh::Edge_index Edge_index;

    typedef std::list<Vertex_index> Curve;
    std::vector<Curve> curves;

    auto collect_constrained_edges = [](const Mesh& mesh, const auto& constraints) {
        std::vector<Edge_index> result;
        for (auto&& edge : mesh.edges()) {
            if (constraints[edge]) {
                //std::cout << "ce:" << static_cast<Point>(mesh.point(mesh.vertex(edge, 0)))
                //    << "->" << static_cast<Point>(mesh.point(mesh.vertex(edge, 1))) << std::endl;
                result.emplace_back(edge);
            }
        }
        return result;
        //std::vector<std::pair<Vertex_index, Vertex_index>> result;
        //for (auto&& edge : mesh.edges()) {
        //    if (constraints[edge]) {
        //        result.emplace_back(
        //            mesh.vertex(edge, 0),
        //            mesh.vertex(edge, 1)
        //        );
        //    }
        //}
        //return result;
    };

    auto constrained_edges = collect_constrained_edges(mesh, constraints);

    auto edge_in_curve = std::map<Edge_index, bool>{};
    for (auto&& edge : constrained_edges) {
        edge_in_curve.emplace(edge, false);
    }
    auto find_free_edge = [&edge_in_curve]() {
        typedef boost::optional<Edge_index> Result;
        auto free_edge = std::find_if(
            begin(edge_in_curve), end(edge_in_curve),
            // find the first pair <key, value> with value = false, i.e. edge is not in curve
            [](auto&& pair) { return !pair.second; }
        );
        if (free_edge == end(edge_in_curve)) return Result{};
        return Result{ free_edge->first };
    };
    auto collect_vertices_until_corner = [&mesh, &edge_in_curve](Vertex_index start, auto out) {
        auto find_single_exit = [&edge_in_curve, &mesh](Vertex_index v) {
            std::cerr << "looking for exit at " << mesh.point(v) << " "
                << std::count_if(begin(edge_in_curve), end(edge_in_curve), [](auto p) { return p.second; })
                << " edges in curves" << std::endl;
            typedef std::pair<Edge_index, Vertex_index> Exit;
            typedef boost::optional<Exit> Result;
            auto result = Result{};
            for (auto&& h : CGAL::halfedges_around_source(v, mesh)) {
                auto p = edge_in_curve.find(mesh.edge(h));
                if (p != end(edge_in_curve)) {
                    if (!p->second) { // edge is not on curve
                        if (result) { // already one exit found
                            return Result{};
                        }
                        else {
                            result = Exit{ p->first, mesh.target(h) };
                        }
                    }
                }
            }
            return result;
        };
        auto next_vertex_found = find_single_exit(start);
        while (next_vertex_found) {
            edge_in_curve[next_vertex_found->first] = true;
            start = next_vertex_found->second;
            *out = start;
            ++out;
            next_vertex_found = find_single_exit(start);
        }
    };
    for (
        auto edge_left = find_free_edge();
        edge_left;
        edge_left = find_free_edge()
        ) {
        auto starting_edge = *edge_left;
        //std::cerr << "#### new curve" << std::endl;
        curves.emplace_back();
        auto& curve = curves.back();
        edge_in_curve[starting_edge] = true;
        curve.emplace_back(mesh.vertex(starting_edge, 0));
        curve.emplace_back(mesh.vertex(starting_edge, 1));
        //std::cerr << "#### new curve backward" << std::endl;
        collect_vertices_until_corner(curve.back(), std::back_inserter(curve));
        //std::cerr << "#### new curve forward" << std::endl;
        collect_vertices_until_corner(curve.front(), std::front_inserter(curve));
    }
    assert(std::all_of(begin(edge_in_curve), end(edge_in_curve), [](auto&& p) {return p.second; }));
    return curves;

}

/** check a curve does goes round a face
(i.e. there is no face that has two edges on the same curve)
otherwise this face should be triangulated */
/**
template <typename Curve>
bool does_curve_round_a_face(const Curve& curve)
{
    auto check_adjacent_faces = [&curve](const Mesh& mesh, auto v1, auto v2) {
        auto other_edges_on_curve = [&curve](const Mesh& mesh, auto h) {
            assert(mesh.is_valid(h));
            if (!mesh.is_border(h)) {
                auto vend = mesh.source(h);
                assert(mesh.point(vend).constraint
                    && mesh.point(vend).constraint.has_curve(&curve));
                assert(mesh.target(h) != vend);
                h = mesh.next(h);
                for (auto v = mesh.target(h); v != vend; h = mesh.next(h)) {
                    if (mesh.point(v).constraint
                        && mesh.point(v).constraint.has_curve(&curve)) return true;
                    v = mesh.target(h);
                }
            }
            return false;
        };
        auto h_on_curve = mesh.halfedge(v1, v2);
        return !(other_edges_on_curve(mesh, h_on_curve) ||
            other_edges_on_curve(mesh, mesh.opposite(h_on_curve)));
    };
    assert(!curve.empty());
    auto pv = begin(curve);
    for (auto qv = next(pv); qv != end(curve); ++qv) {
        for (std::size_t i = 0; i < 2; ++i) {
            auto pi = pv->info[i];
            auto qi = qv->info[i];
            assert(pi.is_on_surface() && qi.is_on_surface());
            assert(pi.surface == qi.surface);
            if (!check_adjacent_faces(*pi.surface, pi.index, qi.index)) return false;
        }
        pv = qv;
    }
    return true;
}


template <typename Vertex_lists>
auto associate_curves(Mesh& mesh1, const Vertex_lists& curves1, Mesh& mesh2, const Vertex_lists& curves2)
{

    typedef typename Vertex_lists::value_type Vertex_list;
    typedef typename Vertex_list::value_type Vertex_index;
    static_assert(std::is_same<typename Mesh::Vertex_index, Vertex_index>::value, "Inconsistent types");

    assert(curves1.size() == curves2.size());
    std::vector<const Vertex_list *> already_associated;
    already_associated.reserve(curves2.size());
    std::map<Vertex_index, Vertex_index> v1tov2, v2tov1;

    for (auto& curve1 : curves1) {
        bool found_corresponding_curve = false;
        for (auto& curve2 : curves2) {
            auto p = &curve2;
            if (std::find(begin(already_associated), end(already_associated), p) == end(already_associated)) {
                if (curve1.size() == curve2.size()) {
                    bool all_points_associated = true;
                    for (auto&& v1 : curve1) {
                        if (v1tov2.count(v1) == 0) {
                            auto P1 = CGAL::Epick::Point_3{ mesh1.point(v1) };
                            auto pv2 = std::find_if(begin(curve2), end(curve2),
                                [&mesh2, &P1](const auto& v2) {
                                return P1 == mesh2.point(v2);
                            });
                            if (pv2 == end(curve2)) {
                                all_points_associated = false;
                                break;
                            }
                            else {
                                assert(mesh1.is_valid(v1));
                                assert(mesh2.is_valid(*pv2));
                                v1tov2[v1] = *pv2;
                            }
                        }
                    }
                    if (all_points_associated) {
                        found_corresponding_curve = true;
                        already_associated.emplace_back(&curve2);
                        break;
                    }
                }
            }
        }
        assert(found_corresponding_curve);
    }

    for (auto&& pair : v1tov2) {
        v2tov1[pair.second] = pair.first;
    }

#ifndef NDEBUG
    auto check_all_points_have_correspondances = [](const auto& curves, const auto& vmap) {
        for (auto&& curve : curves) {
            for (auto&& v : curve) {
                if (vmap.count(v) == 0) return false;
            }
        }
        return true;
    };
    assert(check_all_points_have_correspondances(curves1, v1tov2));
    assert(check_all_points_have_correspondances(curves2, v2tov1));
#endif


    typedef typename Blob::Curve Curve;
    typedef typename Curve::Node Node;
    typedef typename Node::Info Node_info;
    typedef typename Blob::Constraint Constraint;
    typedef typename Constraint::On_curve On_curve;

    std::list<Curve *> curves_on_surface;

    for (auto&& curve1 : curves1) {
        auto pcurve = blob.new_empty_curve();
        assert(pcurve);
        for (auto&& v1 : curve1) {
            const auto v2 = v1tov2[v1];
            assert(mesh1.is_valid(v1));
            assert(mesh2.is_valid(v2));
            pcurve->emplace_back(
                Node_info{ &mesh1, v1 },
                Node_info{ &mesh2, v2 }
            );
            assert(mesh1.point(v1).constraint == mesh2.point(v2).constraint);
            if (!(mesh1.point(v1).constraint)) {
                auto new_constraint = Constraint{
                std::make_shared<On_curve>(pcurve, prev(pcurve->end()))
                };
                mesh1.point(v1).constraint = new_constraint;
                mesh2.point(v2).constraint = new_constraint;
            }
            else {
                assert(mesh1.point(v1).constraint.on_junction());
                assert(false);
            }
        }
        curves_on_surface.emplace_back(pcurve);
    }

    assert(std::all_of(
        begin(curves_on_surface), end(curves_on_surface),
        [](auto p) { return does_curve_round_a_face(*p); }
    ));

    std::vector<Constraint> neighboring_constraints;
    for (auto&& pcurve : curves_on_surface) {
        for (auto&& node : *pcurve) {
            auto v1 = node.info[0].index;
            neighboring_constraints.clear();
            for (auto h1 : CGAL::halfedges_around_source(v1, mesh1)) {
                auto neighbor = mesh1.target(h1);
                if (auto constraint = mesh1.point(neighbor).constraint) {
                    if (auto on_curve = constraint.on_curve()) {
                        if (on_curve->curve != pcurve) {
                            assert(std::find(begin(neighboring_constraints), end(neighboring_constraints), constraint) == end(neighboring_constraints));
                            neighboring_constraints.emplace_back(constraint);
                        }
                    }
                    else {
                        auto on_junction = constraint.on_junction();
                        assert(on_junction);
                    }
                }
            }
            if (neighboring_constraints.size() > 1) {
                for (auto pc = begin(neighboring_constraints); pc != end(neighboring_constraints); ++pc) {
                    for (auto qc = next(pc); qc != end(neighboring_constraints); ++qc) {
                        if (auto weak_link = neighbors_on_same_curve(*pc, *qc)) {
                            assert(weak_link->is_weak());
                            assert(node.constraint());
                            assert(node.constraint().on_curve());
                        }
                    }
                }
            }
        }
    }

}
/**/

