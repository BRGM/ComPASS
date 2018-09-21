#pragma once

template <typename Point_type>
struct Point_id
{
    typedef Surface_type<Point_type> Surface;
    Surface * surface;
    Id_type on_surface_id;
};
