#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

template <typename Mesh>
auto mesh_as_arrays(const Mesh& mesh)
{
    typedef typename Mesh::Point Point;
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

template <typename C3T3, typename Mesh>
void copy_boundary_to_surface_mesh(const C3T3& c3t3,
    typename C3T3::Subdomain_index sd_index,
    Mesh& mesh,
    bool normals_point_outside_of_the_subdomain = true)
{
    typedef typename C3T3::Triangulation Triangulation;
    typedef typename Triangulation::Vertex_handle Vertex_handle;
    typedef typename Mesh::Vertex_index Vertex_index;

    if (normals_point_outside_of_the_subdomain) {
        py::print("Normals point outside subdomain");
    }

    std::map<Vertex_handle, Vertex_index> V;

    int failure_on_first_insertion = 0;
    int failure_on_second_insertion = 0;
    std::array<Vertex_index, 3> indices;
    for (typename C3T3::Facets_in_complex_iterator
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
        fit != end; ++fit)
    {
        typename C3T3::Subdomain_index cell_sd = c3t3.subdomain_index(fit->first);
        typename C3T3::Subdomain_index opp_sd = c3t3.subdomain_index(fit->first->neighbor(fit->second));
        if (cell_sd != sd_index && opp_sd != sd_index) continue;
        if (cell_sd == sd_index && opp_sd == sd_index) continue; // 
        int j = -1;
        for (int i = 0; i < 4; ++i) {
            if (i != fit->second) {
                auto pvi = (*fit).first->vertex(i);
                auto pvi2v = V.find(pvi);
                if (pvi2v == V.end()) {
                    auto v = mesh.add_vertex(pvi->point().point());
                    indices[++j] = v;
                    V[pvi] = v;
                }
                else {
                    indices[++j] = pvi2v->second;
                }
            }
        }
        if (((cell_sd == sd_index) == (fit->second % 2 == 1)) == normals_point_outside_of_the_subdomain)
            std::swap(indices[0], indices[1]);
        auto new_face = mesh.add_face(indices[0], indices[1], indices[2]);
        if (new_face == Mesh::null_face()) {
            ++failure_on_first_insertion; // std::cerr << "First face insertion failed";
            new_face = mesh.add_face(indices[1], indices[0], indices[2]);
            if (new_face == Mesh::null_face()) {
                //py::print(
                //    mesh.halfedge(indices[0], indices[1]) == Mesh::null_halfedge(),
                //    mesh.halfedge(indices[0], indices[2]) == Mesh::null_halfedge(),
                //    mesh.halfedge(indices[1], indices[2]) == Mesh::null_halfedge()
                //);
                ++failure_on_second_insertion;  // std::cerr << " and second too!";
            }
            //std::cerr << std::endl;
        }
    }
    py::print("Failures:", failure_on_first_insertion, failure_on_second_insertion);

}

template <typename Mesh, typename C3T3>
auto collect_all_boundaries_as_surface_meshes(
    const C3T3& c3t3,
    bool normals_point_outside_of_the_subdomain = true)
{
    typedef typename C3T3::Subdomain_index Index;
    std::set<Index> indices;
    for (auto cit = c3t3.cells_in_complex_begin(),
        end = c3t3.cells_in_complex_end();
        cit != end; ++cit) {
        indices.insert(c3t3.subdomain_index(cit));
    }
    auto meshes = py::list{};
    for (auto&& sdi : indices) {
        auto mesh = Mesh{};
        copy_boundary_to_surface_mesh(c3t3, sdi, mesh, normals_point_outside_of_the_subdomain);
        meshes.append(mesh);
    }
    return meshes;
}
