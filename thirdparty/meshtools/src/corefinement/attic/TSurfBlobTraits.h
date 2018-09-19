#pragma once

#include <CGAL/Surface_mesh.h>

#include "Constraint.h"

template <typename Base_point>
struct Constrained_point : Base_point
{
    typedef CGAL::Surface_mesh<Constrained_point> Surface;
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
    Constrained_point(const Base_point& P) :
        Base_point{ P },
        /*id{},*/ constraint{}
    {}
    Constrained_point(Base_point&& P) :
        Base_point{ std::forward<Point>(P) },
        /*id{},*/ constraint{}
    {}
};
