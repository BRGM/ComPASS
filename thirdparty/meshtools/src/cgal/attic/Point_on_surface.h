#pragma once

template <typename Base_point, typename Constraint_type>
struct Constrained_point : Base_point
{
    //typedef Surface_type<Constrained_point> Surface;
    typedef Constraint_type Constraint;
    //Point_id<Point> id;
    Constraint constraint;
    Constrained_point() :
        Base_point{}, /*id{},*/ constraint{}
    {}
    Constrained_point(double x, double y, double z) :
        Base_point{ x, y, z },
        /*id{},*/ constraint{}
    {}
    Constrained_point(const Point& P) :
        Base_point{ P },
        /*id{},*/ constraint{}
    {}
    Constrained_point(Point&& P) :
        Base_point{ std::forward<Point>(P) },
        /*id{},*/ constraint{}
    {}
};
