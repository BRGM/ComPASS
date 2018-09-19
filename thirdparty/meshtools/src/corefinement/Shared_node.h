#pragma once

namespace TSurfBlobTraits
{

    template <typename Surface>
    struct Shared_node
    {
        struct Info
        {
            typedef typename Surface::Vertex_index Vertex_index;
            Surface * surface;
            Vertex_index index;
            bool operator==(const Info& other) const noexcept {
                assert(surface);
                assert(other.surface);
                return surface == other.surface && index == other.index;
            }
            auto point() {
                assert(is_on_surface());
                return surface->point(index);
            }
            auto is_valid() const {
                return surface && (surface->is_valid(index) || index == surface->null_vertex());
            }
            auto is_on_surface() const {
                assert(surface);
                return surface->is_valid(index);
            }
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
        int degree;
        std::array<Info, 3> info;
        Shared_node() = delete;
        Shared_node(Info i1, Info i2) :
            degree{ 2 },
            info{ { i1, i2, Info{} } } {
            assert(i1.surface != i2.surface);
        }
        Shared_node(Info i1, Info i2, Info i3) :
            degree{ 3 },
            info{ { i1, i2, i3 } } {
            assert(
                i1.surface != i2.surface &&
                i1.surface != i3.surface &&
                i2.surface != i3.surface
            );
        }
        Shared_node(Shared_node cn, Info i) :
            degree{ 3 },
            info{ { cn.info[0], cn.info[1], i } } {
            assert(cn.degree == 2);
            assert(
                i.surface != cn.info[0].surface &&
                i.surface != cn.info[1].surface
            );
        }
        //bool is_corner() const { return degree == 3; }
        bool operator==(const Shared_node& other) const noexcept {
            return info == other.info;
        }
        auto point(std::size_t i) noexcept {
            assert(i < 3);
            return info[i].point();
        }
        auto is_valid() const {
            for (int i = 0; i < degree; ++i) {
                if (!info[i].is_valid()) return false;
            }
            return true;
        }
    };

} // namespace TSurfBlobTraits
