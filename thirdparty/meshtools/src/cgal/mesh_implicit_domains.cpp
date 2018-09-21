// This is directly taken from CGAL examples
// cf. https://doc.cgal.org/latest/Mesh_3/examples.html

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh.h>

#include "implicit_functions.h"

using namespace CGAL::parameters;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function> Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Labeled_mesh_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_mesh_domain> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<
    Mesh_domain, CGAL::Default, Concurrency_tag
>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index
> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;
typedef Mesh_criteria::Edge_criteria     Edge_criteria;

#include "mesh-pyutils.h"
#include "mesh_implicit_domains.h"

/** Build the polyline (a circle) which is the intersection
of the sphere and torus function. */
auto sphere_torus_intersection(std::size_t nb_points_on_circle, const double y)
{

    const double r = 5. / 3.;
    const double pi = 2 * acos(0);
    const double dtheta = 2 * pi / (static_cast<double>(nb_points_on_circle));

    typedef typename K::Point_3 Point;
    auto pl = std::vector<Point>{};
    pl.reserve(nb_points_on_circle + 1);
    double theta = 0;
    for (; nb_points_on_circle != 0; --nb_points_on_circle) {
        const double x = r*cos(theta); const double z = r*sin(theta);
        assert(fabs(torus_function(x, y, z)) + fabs(sphere_function<3>(x, y, z))<1E-14);
        if (fabs(torus_function(x, y, z)) + fabs(sphere_function<3>(x, y, z)) > 1E-14)
            std::cerr << "Bad approximation!" << std::endl;
        pl.emplace_back(x, y, z);
        theta += dtheta;
    }
    Point first_point = pl.front(); // pl.emplace_back(pl.front()) fails if pl capacity is saturated
                                    // (e.g. with pl.reserve(nb_points_on_circle))
    pl.emplace_back(first_point); // cycle
    
    return pl;

}

auto sharp_edges(std::size_t nb_points_on_circle = 1000)
{

    const double h = sqrt(2) / 3.;

    std::vector<std::vector<typename K::Point_3>> polylines;
    polylines.emplace_back(
        sphere_torus_intersection(
            nb_points_on_circle, h
        )
    );
    polylines.emplace_back(
        sphere_torus_intersection(
            nb_points_on_circle, -h
        )
    );
    return polylines;

}

py::list mesh_implicit_domains_boundaries()
{
    py::print("Building the model... this may take some time!");
    // Define functions
    Function f1(&torus_function);
    Function f2(&sphere_function<3>);
    Function_vector v;
    v.push_back(f1);
    v.push_back(f2);

    // --- Domain for test case 1
    Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);
    // --- Domain for test case 2
    //std::vector<std::string> vps;
    //vps.push_back("+-");
    //Mesh_domain domain(Function_wrapper(v, vps), K::Sphere_3(CGAL::ORIGIN, 5.*5.));

    //typedef CGAL::Tag_true           Has_features;
    static_assert(std::is_same< CGAL::Tag_true, typename decltype(domain)::Has_features>::value, "Has feature!");

    auto polylines = sharp_edges();
    domain.add_features(begin(polylines), end(polylines));

    // Set mesh criteria
    Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
    Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
    Edge_criteria edge_criteria(0.1); // edge size criteria on sharp edges
    Mesh_criteria criteria{ edge_criteria, facet_criteria, cell_criteria };
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb()
        //, features() // this can be set manually but is true with Mesh_domain_with_polyline_features_3
        //, manifold() // any of the two manifolds options the mesh computation does not end
        //, manifold_with_boundary()
        );

    //py::print("Optimization... this may take even longer!");
    //// Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
    //CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
    //// Exudation
    //CGAL::exude_mesh_3(c3t3, 12);

    py::print("Collect boundaries");
    return collect_all_boundaries_as_surface_meshes<
        CGAL::Surface_mesh<typename CGAL::Epick::Point_3>
        >(c3t3);

}

