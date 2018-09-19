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
        Shared_node(Shared_node cn, Info i) :
            info{ { cn.info[0], cn.info[1], i } } {
            assert(cn.degree() == 2);
            assert(
                i.surface != cn.info[0].surface &&
                i.surface != cn.info[1].surface
            );
        }
        //bool is_corner() const { return degree == 3; }
        bool operator==(const Shared_node& other) const noexcept {
            return info == other.info;
        }
        auto point() noexcept {
            assert(is_consistent());
            return info[0].point();
        }
        auto constraint(std::size_t i) noexcept {
            assert(i<degree());
            return info[i].constraint();
        }
        auto is_valid() const {
            return info[0].is_on_surface() && info[1].is_on_surface() && (
                !info[2].surface || // i.e. degree==2
                info[2].index==info[2].surface->null_vertex() ||
                info[2].is_on_surface()
                );
        }
        auto is_consistent() const {
            assert(is_valid());
            auto P = info[0].base_point();
            return  P == info[1].base_point() && (
                degree() == 2 || is_weak_corner() || P == info[2].base_point()
                );
        }
        auto is_weak_corner() const {
            assert(is_valid());
            return info[2].surface && info[2].index == info[2].surface->null_vertex();
        }
        auto degree() const {
            return info[2].surface ? 3 : 2;
        }
    };

} // namespace TSurfBlobTraits
