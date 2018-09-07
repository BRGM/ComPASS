#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

typedef CGAL::Epick Kernel;

typedef Kernel::Point_3 Point;

//Point id
//(surface, id)
//constraints 
//
//curve_patch = 1 or 2 points if only one point is corner (you don't want to move it)
//connected_point_id, (other_surface, id_on_other_surface)
//
//on deletion -> remove from connected points / test no longer connected


template <typename T>
struct Uid_factory
{
    typedef T Id_type;
    T count;
    Uid_factory() :
        count{ 0 }
    {}
    Uid_factory(const Uid_factory&) = delete;
    Uid_factory& operator=(const Uid_factory&) = delete;
    T make_uid() {
        return count++;
    }
    constexpr T Uid_limit() constexpr const {
        return static_cast<T>(std::numeric_limits<T>::max());
    }
};

// CHECKME: ? Encapsulate id_factory so that we have id{} instead of id{ id_factory.make_uid() }
//Id{
//    static Factory * factory; // = 0
//}

struct Point_with_index: Point
{
    typedef Uid_factory<std::size_t> Factory;
    typedef Factory::Id_type Id;
    static  Factory factory;
    typedef int Degree;
    typedef int Surface_index;
    // CHECKME: ? Encapsulate id_factory so that we have id{} instead of id{ id_factory.make_uid() }
    Id id;
    struct Constraints
    {
        Degree degree; // 0 initialisation
        std::array<Surface_index, 3> surfaces;
        auto begin() { return surfaces.data(); }
        auto end() { return surfaces.data() + degree; }
        bool has_constraint(Surface_index si) {
            return std::find(begin(), end(), si) != end();
        }
        void add_constraint(Surface_index si) {
            assert(!has_constraint(si));
            assert(degree < surfaces.size());
            surfaces[degree] = si;
            ++degree;
            if (degree > 1) {
                std::sort(begin(), end());
            }
        }
    };
    Constraints constraints;
    Point_with_index() :
        Point{},
        id{ factory.make_uid() }
    {}
    Point_with_index(double x, double y, double z) :
        Point{ x, y, z },
        id{ factory.make_uid() }
    {}
    Point_with_index(const Point& P) :
        Point{ P },
        id{ factory.make_uid() }
    {}
    Point_with_index(Point&& P) :
        Point{ std::forward<Point>(P) },
        id{ factory.make_uid() }
    {}
};

Point_with_index::Factory Point_with_index::factory;

typedef CGAL::Surface_mesh<Point_with_index> Mesh;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

auto mesh_as_arrays(const Mesh& mesh)
{
    typedef std::size_t New_index;
    std::vector<New_index> reindex;
    reindex.resize(mesh.num_vertices());
    auto nv = New_index{ 0 };
    for (auto&& v : mesh.vertices()) {
        reindex[v] = nv++;
    }
    auto vertices = py::array_t<double, py::array::c_style>{ 
        { static_cast<std::size_t>(nv), static_cast<std::size_t>(Point::Ambient_dimension::value) }
    };
    {
        auto p = reinterpret_cast<Point*>(vertices.request().ptr);
        static_assert(sizeof(Point) == Point::Ambient_dimension::value * sizeof(double), "Inconsistent sizes in memory!");
        for (auto&& v : mesh.vertices()) {
            *(p++) = mesh.point(v);
        }
    }
    auto triangles = py::array_t<New_index, py::array::c_style>{ 
        { static_cast<std::size_t>(mesh.number_of_faces()), static_cast<std::size_t>(3) }
    };
    {
        auto p = reinterpret_cast<New_index*>(triangles.request().ptr);
        for (auto&& f : mesh.faces()) {
            assert(mesh.degree(f) == 3);
            for (auto&& v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
                *(p++) = reindex[v];
            }
        }
    }
    return py::make_tuple(vertices, triangles);
}

auto test()
{
    Mesh tm1;
    Mesh::Vertex_index u = tm1.add_vertex(Point(0, 1, 0));
    Mesh::Vertex_index v = tm1.add_vertex(Point(0, 0, 0.2));
    Mesh::Vertex_index w = tm1.add_vertex(Point(1, 0, 0));
    tm1.add_face(u, v, w);
    Mesh tm2;
    u = tm2.add_vertex(Point(  0, 0.4, -1));
    v = tm2.add_vertex(Point(  1, 0.5, -1));
    w = tm2.add_vertex(Point(0.5, 0.6,  1));
    tm2.add_face(u, v, w);
    CGAL::Polygon_mesh_processing::corefine(tm1, tm2, true);
    std::ofstream os("test.off");
    CGAL::write_off(os, tm1);
    os.close();
    return py::make_tuple(mesh_as_arrays(tm1), mesh_as_arrays(tm2));
}

PYBIND11_MODULE(Corefinement, module)
{

    module.doc() = "pybind11 with CGAL surface corefinement (quick and dirty!!!)";

    module.def("test", &test);

}
