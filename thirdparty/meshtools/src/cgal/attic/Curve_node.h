#pragma once

template <typename Surface>
struct Constrained_node
{
    struct Info
    {
        typedef typename Surface::Vertex_index Vertex_index;
        Surface * surface;
        Vertex_index index;
        bool operator==(const Surface_point_index& other) const noexcept {
            assert(surface);
            assert(other.surface);
            return surface == other.surface && index == other.index;
        }
        auto point() {
            assert(surface);
            assert(surface->is_valid(index));
            return surface->point(index);
        }
    };
    int degree;
    std::array<Info, 3> info;
    Curve_node() = delete;
    Curve_node(Info i1, Info i2) :
        degree{ 2 },
        info{ { i1, i2, Info{} } } {
        assert(i1.surface != i2.surface);
    }
    Curve_node(Info i1, Info i2, Info i3) :
        degree{ 3 },
        info{ { i1, i2, i3 } } {
        assert(
            i1.surface != i2.surface &&
            i1.surface != i3.surface &&
            i2.surface != i3.surface
        );
    }
    Curve_node(Curve_node cn, Info i) :
        degree{ 3 },
        info{ { cn.info[0], cn.info[1], i } } {
        assert(cn.degree == 2);
        assert(
            i.surface != cn.info[0].surface &&
            i.surface != cn.info[1].surface
        );
    }
    //bool is_corner() const { return degree == 3; }
    bool operator==(const Curve_node& other) const noexcept {
        return info == other.info;
    }
    auto point(std::size_t i) noexcept {
        assert(i < 3);
        return info[i].point();
    }
};
