#pragma once

#include <CGAL/Surface_mesh.h>

template <typename BlobPoint>
struct TSurfBlob
{

    typedef BlobPoint Point;
    typedef CGAL::Surface_mesh<Point> TSurf;
    typedef std::list<Point> Curve;
    
    typedef std::pair<TSurf *, TSurf *> Surface_intersection;

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
    std::multimap<Surface_intersection, Curve *> intersections;

    auto new_empty_curve() {
        return curves_factory.new_curve();
    }

    auto new_intersection(TSurf * S1, TSurf * S2) {
        auto curve = new_empty_curve();
        assert(curve);
        intersections.insert({ { S1, S2 }, curve });
        return curve;
    }

};
