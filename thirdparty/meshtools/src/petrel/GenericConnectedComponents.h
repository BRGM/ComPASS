#pragma once

#include <cassert>
#include <set>
#include <list>
#include <map>

template <typename ElementHandle>
class ConnectedComponents
{
public:
    typedef ElementHandle Element_handle;
    typedef std::vector<Element_handle> Component;
protected:
    typedef std::vector<Component> Components;
    Components components;
public:
    typedef typename Components::const_iterator Component_iterator;
    template <typename InputIterator, typename NeighborhoodFactory>
    ConnectedComponents(InputIterator first, InputIterator last, const NeighborhoodFactory& neighborhood_factory) :
        components{}
    {
        typedef std::size_t Component_id;
        typedef std::map<Element_handle, Component_id> Components_map;
        std::vector<Component_id> connections;
        Components_map components_map;
        Component discovered_neighbors;
        std::vector<typename Components_map::iterator> existing_components;
        // FIXME: shoud have a lock mechanism for parrallel processing
        for (auto p = first; p != last; ++p) {
            discovered_neighbors.clear();
            existing_components.clear();
            auto handle = Element_handle{ p };
            auto element_component = components_map.find(handle);
            if (element_component == components_map.end()) {
                auto result = components_map.insert({ handle, Component_id{ connections.size() } });
                assert(result.second);
                element_component = result.first;
            }
            auto neighborhood = neighborhood_factory(handle);
            for (int i = 0; i < neighborhood.size(); ++i) {
                if (neighborhood.is_connected(i)) {
                    auto neighbor = neighborhood(i);
                    auto neighbor_component = components_map.find(neighbor);
                    if (neighbor_component != components_map.end())
                        existing_components.push_back(neighbor_component);
                    else
                        discovered_neighbors.push_back(neighbor);
                }
            }
            auto component = element_component->second;
            if (existing_components.empty()) {
                connections.emplace_back(component);
            }
            else {
                existing_components.push_back(element_component);
                for (auto&& other_component : existing_components) {
                    component = std::min(component, connections[other_component->second]);
                }
                assert(component < connections.size());
                element_component->second = component;
                for (auto&& other_component : existing_components) {
                    connections[other_component->second] = component;
                }
            }
            for (auto&& neighbor : discovered_neighbors) {
                components_map.insert({neighbor, component});
            }
        }
        std::vector<Component_id> path;
        auto next_component = [&connections](Component_id ci) { return connections[ci]; };
        for (auto p = crbegin(connections); p != crend(connections); ++p) {
            auto c = connections[*p];
            if (c != next_component(c)) {
                path.clear();
                path.push_back(*p);
                for (; c != next_component(c); c = next_component(c)) {
                    path.push_back(c);
                }
                auto target = path.back();
                assert(target == next_component(target));
                for (auto c_on_path : path) {
                    connections[c_on_path] = target;
                }
            }
        }
        std::set<Component_id> components_ids{ connections.cbegin(), connections.cend() };
        components = Components{ components_ids.size() };
        std::map<std::size_t, std::size_t> component_final_id;
        std::size_t pos = 0;
        for (auto id : components_ids) {
            component_final_id[id] = pos;
            ++pos;
        }
        for (auto p = first; p != last; ++p) {
            auto handle = Element_handle{ p };
            components[component_final_id[connections[components_map[handle]]]].push_back(handle);
        }
        for (auto&& component : components) {
            component.shrink_to_fit();
        }
    }
    auto nb_components() const { return components.size(); }
    Component_iterator begin() const { return components.begin(); }
    Component_iterator end() const { return components.end(); }
};

//template <typename ElementHandle, typename InputIterator,
//    typename ValidityTest, typename ConnexionTest,
//    typename Neighborhood>
//    ConnectedComponents<ElementHandle>::ConnectedComponents(
//        InputIterator first, InputIterator last,
//        const Neighborhood& neighborhood,
//        const ValidityTest& is_valid,
//        const ConnexionTest& are_connected
//    )
////
//template <typename Complex >
//auto connected_components(const Complex& complex)
//{
//    return ConnectedComponents<Complex>{ complex };
//}
//
//template <typename Complex, typename CellId = int >
//auto convert_components(const Complex& complex, const ConnectedComponents<Complex>& components)
//{
//    typedef typename Complex::Cell_handle Cell_handle;
//    std::map<Cell_handle, CellId> cmap;
//    auto ci = CellId{ 0 };
//    for (auto cell = complex.cells_in_complex_begin(); cell != complex.cells_in_complex_end(); ++cell) {
//        Cell_handle handle = cell;
//        cmap[handle] = ci;
//        ++ci;
//    }
//    std::list< std::vector< CellId > > converted_components;
//    for (auto&& component : components) {
//        converted_components.emplace_back();
//        auto& new_component = converted_components.back();
//        new_component.reserve(component.size());
//        for (auto&& cell : component) {
//            assert(complex.is_in_complex(cell));
//            new_component.emplace_back(cmap[cell]);
//        }
//    }
//    return converted_components;
//}