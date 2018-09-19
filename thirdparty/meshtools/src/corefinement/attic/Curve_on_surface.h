#pragma once

template <typename Node_type>
struct Curve : std::list<Node_type>
{
    typedef iterator Position;
    typedef Node_type Node;
    Curve() = delete;
    Curve(const Curve&) = delete;
    Curve& operator=(const Curve&) = delete;
    // This is just an empty struct to restrict the creation of curves object
    struct Creation_tag {};
    Curve(const Creation_tag) :
        std::list<Curve_node<Surface>>{}
    {}
};

