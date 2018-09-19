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
    };

    Curves_factory curves_factory;

    auto new_empty_curve() {
        return curves_factory.new_curve();
    }

    auto split(On_curve& on_curve)
    {
        typedef typename Constraint::On_corner::Link Link;
        assert(on_curve.is_valid());
        auto head = on_curve.curve;
        auto tail = new_empty_curve();
        head->splice(*tail, end(*tail), on_curve.position, end(*head));
        assert(tail->size() >= 1);
        head.emplace_back(tail.front());
        auto p = begin(*tail);
        while (next(p) != end(*tail)) {
            assert(p->point(0).constraint);
            auto ocp = p->point(0).constraint.on_curve();
            assert(ocp);
            assert(ocp->curve == head);
            ocp->curve = tail;
            ocp->position = p;
            assert(p->point(1).constraint.on_curve());
            assert(p->point(1).constraint.on_curve()->curve == tail);
            assert(p->point(1).constraint.on_curve()->position == p);
            assert(p->degree() == 2);
            ++p;
        }
        auto latest_constraint = p->point(1).constraint;
        if (auto ocp = latest_constraint.on_corner()) {
            for (auto& link : ocp->links) {
                assert(link.is_valid());
                assert(!link.is_weak());
                if (link.first.curve == head) {
                    link.first.curve = tail;
                    link.first.position = p;
                }
                assert(link.second.curve != head);
                assert(link.is_valid());
            }
        }
        else {
            auto ocup = latest_constraint.on_curve();
            assert(ocup);
            assert(ocup->is_on_extremity());
            ocup->curve = tail;
            ocup->position = p;
            assert(p->point(1).constraint.on_curve());
            assert(p->point(1).constraint.on_curve()->curve == tail);
            assert(p->point(1).constraint.on_curve()->position == p);
        }
        assert(p->point(1).constraint.on_corner());
            //if (p)
            //    assert(head->size() >= 1);
        // original constraint is deactivated
        on_curve.curve = nullptr;
        auto link = Link{
            { head, prev(head->end()) },
            { tail, tail->begin() }
        };
        assert(link.is_valid());
        return link;
    }

};
