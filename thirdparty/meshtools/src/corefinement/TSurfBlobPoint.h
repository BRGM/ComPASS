#pragma once

#include <CGAL/Surface_mesh.h>

#include "Shared_node.h"
#include "Curve.h"
#include "Constraint.h"

namespace BlobTraits = TSurfBlobTraits;

template <typename BasePoint>
struct TSurfBlobPoint : BasePoint
{
    typedef BasePoint Base_point;
    typedef CGAL::Surface_mesh<TSurfBlobPoint> TSurf;
    typedef BlobTraits::Shared_node<TSurf> Shared_node;
    typedef BlobTraits::Curve<Shared_node> Curve;
    typedef BlobTraits::Constraint<Curve> Constraint;
    //Point_id<Point> id;
    Constraint constraint;
    TSurfBlobPoint() :
        Base_point{}, /*id{},*/ constraint{}
    {}
    TSurfBlobPoint(double x, double y, double z) :
        Base_point{ x, y, z },
        /*id{},*/ constraint{}
    {}
    TSurfBlobPoint(const Base_point& P) :
        Base_point{ P },
        /*id{},*/ constraint{}
    {}
    TSurfBlobPoint(Base_point&& P) :
        Base_point{ std::forward<Base_point>(P) },
        /*id{},*/ constraint{}
    {}
};

