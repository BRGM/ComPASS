#pragma once

namespace TSurfBlobTraits
{

    template <typename Surface>
    struct Shared_node
    {
        struct Info
        {
            typedef typename Surface::Point Surface_point;
            typedef typename Surface_point::Base_point Surface_base_point;
            typedef typename Surface_point::Constraint Surface_constraint;
            typedef typename Surface::Vertex_index Vertex_index;
            Surface * surface;
            Vertex_index index;
            bool operator==(const Info& other) const noexcept {
                assert(surface);
                assert(other.surface);
                return surface == other.surface && index == other.index;
            }
            Surface_point& point() {
                assert(is_on_surface());
                return surface->point(index);
            }
            const Surface_point& point() const {
                assert(is_on_surface());
                return surface->point(index);
            }
            Surface_base_point& base_point() {
                return point();
            }
            const Surface_base_point& base_point() const {
                return point();
            }
            Surface_constraint& constraint() {
                return point().constraint;
            }
            const Surface_constraint& constraint() const {
                return point().constraint;
            }
            auto is_valid() const {
                return surface && (surface->is_valid(index) || index == surface->null_vertex());
            }
            auto is_on_surface() const {
                assert(surface);
                return surface->is_valid(index);
            }
            Info() :
                surface(nullptr),
                index{} 
            {}
            Info(Surface * S) :
                surface(S),
                index{ surface->null_vertex() } {
                assert(!is_on_surface());
            }
            Info(Surface * S, Vertex_index v) :
                surface(S),
                index{ v } {
                assert(is_on_surface());
            }
        };
        std::array<Info, 3> info;
        Shared_node() = delete;
        Shared_node(Info i1, Info i2) :
            info{ { i1, i2, Info{} } } {
            assert(i1.surface != i2.surface);
        }
        Shared_node(Info i1, Info i2, Info i3) :
            info{ { i1, i2, i3 } } {
            assert(
                i1.surface != i2.surface &&
                i1.surface != i3.surface &&
                i2.surface != i3.surface
            );
        }
        Shared_node(Shared_node node, Info i) :
            info{ { node.info[0], node.info[1], i } } {
            assert(node.degree() == 2);
            assert(
                i.surface != node.info[0].surface &&
                i.surface != node.info[1].surface
            );
        }
        //bool is_corner() const { return degree == 3; }
        bool operator==(const Shared_node& other) const noexcept {
            return info == other.info;
        }
        auto base_point() noexcept {
            assert(is_geometrically_consistent());
            if(info[0].is_on_surface()) return info[0].base_point();
            return info[1].base_point();
        }
        auto point() noexcept {
            assert(is_consistent());
            if (info[0].is_on_surface()) return info[0].point();
            return info[1].point();
        }
        auto constraint() noexcept {
            return point().constraint;
        }
        auto is_valid() const {
            if(!info[0].is_valid()) return false;
            if(!info[1].is_valid()) return false;
            if (info[2].surface) {
                if (!info[2].is_valid()) return false;
            }
            if (info[0].surface == info[1].surface) return false;
            if (info[0].surface == info[2].surface) return false;
            if (info[1].surface == info[2].surface) return false;
            if (!info[0].is_on_surface()) {
                if (degree() != 3) return false;
                if (!(info[1].is_on_surface() && info[2].is_on_surface())) return false;
            }
            if (!info[1].is_on_surface()) {
                if (degree() != 3) return false;
                if (!(info[0].is_on_surface() && info[2].is_on_surface())) return false;
            }
            if (degree()==3 && !info[2].is_on_surface()) {
                if (!(info[0].is_on_surface() && info[1].is_on_surface())) return false;
            }
            return true;
        }
        auto is_consistent() const {
            assert(is_valid());
            if (degree() == 2) return info[0].point() == info[1].point();
            std::set<typename Info::Surface_point> pts;
            for (auto&& i : info) {
                if (i.is_on_surface()) pts.insert(i.point());
            }
            return pts.size()==1;
        }
        auto is_geometrically_consistent() const {
            assert(is_valid());
            if (degree() == 2) return info[0].base_point() == info[1].base_point();
            std::set<typename Info::Surface_base_point> pts;
            for (auto&& i : info) {
                if (i.is_on_surface()) pts.insert(i.base_point());
            }
            return pts.size() == 1;
        }
        auto is_weak_corner() const {
            assert(is_valid());
            return info[2].surface && (
                info[0].index == info[0].surface->null_vertex() ||
                info[1].index == info[1].surface->null_vertex() ||
                info[2].index == info[2].surface->null_vertex()
            );
        }
        auto degree() const {
            return info[2].surface ? 3 : 2;
        }
    };

} // namespace TSurfBlobTraits
