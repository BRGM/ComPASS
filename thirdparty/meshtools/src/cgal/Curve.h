#pragma once

namespace TSurfBlobTraits
{

    template <typename Node_type>
    struct Curve : std::list<Node_type>
    {
        typedef typename std::list<Node_type>::iterator Position;
        typedef Node_type Node;
        Curve() = default;
        Curve(const Curve&) = delete;
        Curve& operator=(const Curve&) = delete;
        // This is just an empty struct to restrict the creation of curves object
        bool is_valid() const {
            if(this->size() < 2) return false;
            auto S1 = this->front().info[0].surface;
            auto S2 = this->front().info[1].surface;
            for (auto& node : *this) {
                if (!node.is_valid()) return false;
                if (node.info[0].surface != S1) return false;
                if (node.info[1].surface != S2) return false;
            }
            for (auto pnode = next(this->begin()); pnode != prev(this->end()); ++pnode) {
                if (pnode->degree() != 2) return false;
            }
            return true;
        }
        auto incident_surfaces() const {
            assert(is_valid());
            return std::make_pair(
                this->front().info[0].surface,
                this->front().info[1].surface
            );
        }
        template <typename Surface>
        auto has_incident_surface(Surface * S) const {
            assert(is_valid());
            return (
                this->front().info[0].surface == S ||
                this->front().info[1].surface == S;
            );
        }
        auto share_surface(const Curve& other) const {
            auto this_surfaces = incident_surfaces();
            auto other_surfaces = other.incident_surfaces();
            if (this_surfaces.first == other_surfaces.first) return this_surfaces.first;
            if (this_surfaces.first == other_surfaces.second) return this_surfaces.first;
            if (this_surfaces.second == other_surfaces.first) return this_surfaces.second;
            if (this_surfaces.second == other_surfaces.second) return this_surfaces.second;
            return nullptr;
        }
        auto complementary_surface(const Curve& other) const {
            auto common_surface = share_surface(other);
            assert(common_surface);
            auto other_surfaces = other.incident_surfaces();
            if (other_surfaces.first == common_surface) return other_surfaces.second;
            assert(other_surfaces.second == common_surface);
            return other_surfaces.first;
        }
    };

} // namespace TSurfBlobTraits
