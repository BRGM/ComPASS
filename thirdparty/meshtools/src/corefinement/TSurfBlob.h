#pragma once

template <typename BlobPoint>
struct TSurfBlob
{

    typedef BlobPoint Point;
    typedef typename Point::TSurf TSurf;
    typedef typename Point::Curve Curve;
    typedef typename Point::Constraint Constraint;
    typedef typename Constraint::On_curve On_curve;

    struct Curves_factory
    {
    private:
        std::vector<std::shared_ptr<Curve>> pointers;
    public:
        Curve * new_curve() {
            pointers.emplace_back(
                std::make_shared<Curve>()
            );
            return pointers.back().get();
        }
        auto already_created() const noexcept {
            return pointers.size();
        }
    };

    Curves_factory curves_factory;

    auto new_empty_curve() {
        return curves_factory.new_curve();
    }

    void reassociate(Curve& curve, Curve * old)
    {
        for (auto p = begin(curve); p != end(curve); ++p) {
            assert(p->is_geometrically_consistent());
            if (auto on_curve = p->constraint(0).on_curve()) {
                assert(on_curve->curve == old);
                on_curve->curve = &curve;
                on_curve->position = p;
            }
            else {
                auto on_junction = p->constraint(0).on_junction();
                assert(on_junction);
                assert(on_junction->is_valid());
                assert(next(p) == end(*old));
                assert(on_junction->has_incoming_curve(old));
                for (auto&& junction : on_junction->junctions) {
                    if (junction.first.curve == old) {
                        junction.first.curve = &curve;
                        junction.first.position = p;
                    }
                }
            }
        }
    }

    auto split(On_curve& on_curve, TSurf * S)
    {
        assert(on_curve.is_valid());
        assert(S);
        assert(!on_curve.is_on_extremity());
        auto head = on_curve.curve;
        assert(head->is_valid());
        auto tail = new_empty_curve();
        head->splice(*tail, end(*tail), next(on_curve.position), end(*head));
        assert(on_curve.position == prev(head->end()));
        assert(on_curve.is_on_extremity());
        assert(on_curve->position->degee() == 2);
        on_curve->position->info[2] = S;
        assert(head->back().is_weak_corner());
        tail->emplace_front(head.back());
        reassociate(*tail, head);
        assert(head->is_valid());
        assert(tail->is_valid());
        return tail;
    }

    auto split(
        Curve * head, typename Curve::Position p, 
        TSurf * S0, typename TSurf::Vertex_index v0,
        TSurf * S2, typename TSurf::Vertex_index v2
    ) {
        assert(head->is_valid());
        assert(next(p)!=end(*head));
        assert((p->contraint(0).on_curve() && p->contraint(0).on_curve()->curve == head) ||
            (p->contraint(0).on_junction() && p->contraint(0).on_junction()->has_outgoing_curve(head))
        );
        assert((next(p)->contraint(0).on_curve() && next(p)->contraint(0).on_curve()->curve == head) ||
            (next(p)->contraint(0).on_junction() && next(p)->contraint(0).on_junction()->has_incoming_curve(head))
        );
        assert(S0);
        assert(S0->is_valid(v0));
        assert(S2);
        assert(S2->is_valid(v2));
        assert(S0 != S2);
        typedef typename Curve::Node Curve_node;
        typedef typename Curve_node::Info Node_info;
        auto make_new_node = [p, S0, v0, S2, v2]() {
            if (p->info[0].surface == S0) {
                return Curve_node{
                    Node_info{ S0, v0 },
                    Node_info{ p->info[1].surface },
                    Node_info{ S2, v2 }
                };
            }
            assert(p->info[1].surface == S0);
            return Curve_node{
                Node_info{ p->info[0].surface },
                Node_info{ S0, v0 },
                Node_info{ S2, v2 }
            };
        };
        auto new_node = make_new_node();
        assert(new_node.is_on_weak_corner());
        auto tail = new_empty_curve();
        head->splice(*tail, end(*tail), next(p), end(*head));
        head->emplace_back(new_node);
        tail->emplace_front(new_node);
        reassociate(*tail, head);
        assert(head->is_valid());
        assert(tail->is_valid());
        return tail;
    }

};
