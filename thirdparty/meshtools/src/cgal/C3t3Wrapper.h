#pragma once

#include <iostream>
#include <fstream>
#include <utility>
#include <iterator>
#include <memory>
#include <functional>
#include <sstream>

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Mesh_domain_with_polyline_features_3.h"
#include "CGAL/Mesh_triangulation_3.h"
#include "CGAL/Mesh_complex_3_in_triangulation_3.h"
#include <CGAL/IO/File_binary_mesh_3.h>

#include "DummyDomain.h"
#include "BufferInfo.h"
//#include "DistanceTrees.h"

typedef CGAL::Epick Kernel;
typedef int SubdomainIndex;
typedef std::pair<SubdomainIndex, SubdomainIndex> SurfacePatchIndex;

typedef ExperimentalOracle::DummyDomain<Kernel, SubdomainIndex, SurfacePatchIndex> Domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Domain> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain, Mesh_domain::R, CGAL::Sequential_tag>::type Triangulation;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Triangulation> C3t3;
//typedef ExperimentalBRGM::DistanceTrees<C3t3> DTree;
//typedef DTree::Point DTP;

// Argument Dependent Lookup for operator<<
namespace CGAL
{

	std::istream& operator >> (std::istream& is, SurfacePatchIndex& spi)
	{

		for (char c = ' '; c != '('; is >> c) {}
		is >> spi.first;
		for (char c = ' '; c != ','; is >> c) {}
		is >> spi.second;
		for (char c = ' '; c != ')'; is >> c) {}
		return is;

	}

}

bool load_c3t3(const std::string& filename, C3t3& c3t3)
{

	std::cout << "Extracting data from file: " << filename << std::endl;
	//std::ifstream file(argv[1], std::ios_base::binary | std::ios_base::in);
	//CGAL::set_binary_mode(file);
	std::ifstream file(filename, std::ios_base::in);
	if (file.good()) {
		CGAL::set_ascii_mode(file);
		file >> c3t3;
		std::cout << "Complex statistics:" << std::endl;
		std::cout << c3t3.number_of_cells() << " cells" << std::endl;
		std::cout << c3t3.number_of_facets() << " facets" << std::endl;
		file.close();
		// output to medit
		std::ofstream medit_file("out.mesh");
		c3t3.output_to_medit(medit_file);
		medit_file.close();
	}
	else {
		file.close();
		std::cerr << "Something went wrong loading c3t3." << std::endl;
        throw std::runtime_error("Could not extract C3T3 from " + filename);
        return false;
	}
	file.close();
	return true;

}

bool binary_load_c3t3(const std::string& filename, C3t3& c3t3)
{

    std::cout << "Extracting data from file: " << filename << std::endl;
    std::ifstream file(filename, std::ios_base::binary | std::ios_base::in);
    bool ok = file.good();
    if (ok) {
        ok = CGAL::Mesh_3::load_binary_file(file, c3t3);
        std::cout << "Complex statistics:" << std::endl;
        std::cout << c3t3.number_of_cells() << " cells" << std::endl;
        std::cout << c3t3.number_of_facets() << " facets" << std::endl;
    }
    file.close();
    if(!ok) {
        std::cerr << "Something went wrong loading c3t3." << std::endl;
        throw std::runtime_error("Could not extract C3T3 from " + filename);
        return false;
    }
    return true;

}

struct Pt
{
	double x;
	double y;
	double z;
	auto as_string() const {
		std::ostringstream os;
		os << "(" << x << "," << y << "," << z << ")";
		return os.str();
	}
};

template<typename R, typename I>
auto recast_function(R(*f)(I)) {
	return std::function<R(I)>{f};
}

template <typename T>
class Id_generator
{
protected:
	T id;
public:
	Id_generator(T start = 0) :
		id{ start } {}
	Id_generator(const Id_generator&) = delete;
	Id_generator& operator=(const Id_generator&) = delete;
	T operator()() { return id++; }
	T status() const { return id; }
};

class C3t3Wrapper
{
public:
	constexpr static std::size_t dim = 3;
	typedef double Field_type;
	typedef int Id_type;
	typedef std::array<Field_type, dim> Vertex;
	typedef std::array<Id_type, dim + 1> Simplex;
	typedef std::array<Id_type, dim> Facet;
	typedef Id_type Domain_tag;
	typedef std::array<Id_type, 2> Facet_tag;
protected:
	C3t3 c3t3;
	std::vector<Vertex> vertices;
	std::vector<Simplex> cells;
	std::vector<Facet> facets;
	std::vector<Domain_tag> domain_tags;
	std::vector<Facet_tag> facet_tags;
	//std::shared_ptr<DTree> dtree;
	template <typename Container>
	void reset_container(Container& container, std::size_t target_size) {
		container.clear();
		container.reserve(target_size);
		assert(container.empty());
	}
	void reset_all_containers() {
		const C3t3& complex = c3t3;
		const Triangulation& triangulation = complex.triangulation();
		reset_container(vertices, triangulation.number_of_vertices());
		const std::size_t nbcells = complex.number_of_cells_in_complex();
		reset_container(cells, nbcells);
		reset_container(domain_tags, nbcells);
		const std::size_t nbfacets = complex.number_of_facets_in_complex();
		reset_container(facets, nbfacets);
		reset_container(facet_tags, nbfacets);
	}
	void update_info() {
		reset_all_containers();
		const C3t3& complex = c3t3;
		const Triangulation& triangulation = complex.triangulation();
		typedef std::map<typename Triangulation::Vertex_handle, Id_type> Vertices_map;
		Vertices_map vmap;
		{
			Id_generator<Id_type> genid;
			for (auto p = triangulation.finite_vertices_begin(); p != triangulation.finite_vertices_end(); ++p) {
				vmap[p] = genid();
			}
			assert(genid.status() == triangulation.number_of_vertices());
		}
		typedef Triangulation::Vertex_iterator::value_type Triangulation_vertex;
		std::transform(
			triangulation.finite_vertices_begin(), triangulation.finite_vertices_end(),
			std::back_inserter(vertices),
			[](const Triangulation_vertex& p) -> Vertex {
			auto P = p.point();
			return Vertex{ { P.x(), P.y(), P.z() } };
		});
		//std::vector<Field_type>& tmp = vertices_tmp;
		//std::for_each(
		//	triangulation.finite_vertices_begin(), triangulation.finite_vertices_end(),
		//	[&tmp](const Triangulation_vertex& p) {
		//	auto P = p.point();
		//	std::copy(P.cartesian_begin(), P.cartesian_end(), std::back_inserter(tmp));
		//});
		typedef C3t3::Cells_in_complex_iterator::value_type Complex_cell;
		std::transform(
			complex.cells_in_complex_begin(), complex.cells_in_complex_end(),
			std::back_inserter(cells),
			[&vmap](const Complex_cell& cell) -> Simplex {
			return Simplex{ {
					vmap.at(cell.vertex(0)),
					vmap.at(cell.vertex(1)),
					vmap.at(cell.vertex(2)),
					vmap.at(cell.vertex(3))
				} };
		});
		std::transform(
			complex.cells_in_complex_begin(), complex.cells_in_complex_end(),
			std::back_inserter(domain_tags),
			[](const Complex_cell& cell) -> Domain_tag { 
				return cell.subdomain_index(); 
		});
		typedef C3t3::Facets_in_complex_iterator::value_type Complex_facet;
		std::transform(
			complex.facets_in_complex_begin(), complex.facets_in_complex_end(),
			std::back_inserter(facets),
			[&vmap](const Complex_facet& facet) -> Facet {
			const auto cell = facet.first;
			auto nodes = Facet{};
			int pos = 0;
			for (int i = 0; i < 4; ++i) {
				if (i != facet.second) {
					nodes[pos] = vmap.at(cell->vertex(i));
					++pos;
				}
			}
			return nodes;
		});
		std::transform(
			complex.facets_in_complex_begin(), complex.facets_in_complex_end(),
			std::back_inserter(facet_tags),
			[complex](const Complex_facet& facet) -> Facet_tag {
			const auto patch_index = complex.surface_patch_index(facet);
			return Facet_tag{ { patch_index.first, patch_index.second } };
		});
		assert(vertices.size() == triangulation.number_of_vertices());
		assert(cells.size() == complex.number_of_cells());
		assert(facets.size() == complex.number_of_facets());
		assert(facet_tags.size() == complex.number_of_facets());
		//dtree = std::make_shared<DTree>(complex);
	}
public:
    void binary_reload(const std::string& filename) {
        binary_load_c3t3(filename, c3t3);
        update_info();
    }
    void reload(const std::string& filename) {
        load_c3t3(filename, c3t3);
        update_info();
    }
    //   auto squared_distance_functor(Id_type f1, Id_type f2) const {
	//	auto tree = std::shared_ptr<DTree>{ dtree };
	//	auto spi = SurfacePatchIndex{ f1, f2 };
	//	return [tree, spi](const Pt& P) -> double {
	//		static_assert(sizeof(DTree::Point) == sizeof(Pt), "Incompatible Point types !");
	//		return tree->squared_distance(reinterpret_cast<const DTree::Point&>(P), spi);
	//	};
	//}
    auto vertices_buffer() const { return BufferInfo<double, 3>{ vertices }; }
    auto cells_buffer() const { return BufferInfo<Id_type, dim + 1>{ cells }; }
    auto domain_tags_buffer() const { return BufferInfo<Domain_tag>{ domain_tags }; }
    auto facets_buffer() const { return BufferInfo<Id_type, dim>{ facets }; }
    auto facet_tags_buffer() const { return BufferInfo<Id_type, 2>{ facet_tags }; }
    const C3t3& get() const {
        return c3t3;
    }
};

//template<typename Result, typename Input>
//auto apply(std::function<Result(Input)> F, py::array_t<double, py::array::c_style> a)
//{
//    typedef typename std::remove_reference<Input>::type Input_type;
//    typedef py::array_t<Result, py::array::c_style> Result_array;
//    if (a.ndim() == 0) return Result_array{ {} };
//    std::vector<std::size_t> result_shape;
//    // IMPROVE: set as boolean template
//    if (a.itemsize() == sizeof(Input_type)) {
//        auto shape = a.request().shape;
//        result_shape.reserve(shape.size());
//        for (auto&& n : shape) {
//            result_shape.push_back(n);
//        }
//    }
//    else if (a.shape(a.ndim() - 1) * a.itemsize() == sizeof(Input_type)) {
//        const auto info = a.request();
//        result_shape = std::vector<std::size_t>{ info.shape.begin(), info.shape.begin() + a.ndim() - 1 };
//    }
//    else {
//        throw std::runtime_error("Incompatible size in memory!");
//    }
//    auto result = Result_array{ result_shape };
//    // WARNING: If Input_type is not const cast will fail as a.data is const: use a.mutable_data
//    auto begin_input = reinterpret_cast<Input_type *>(a.data(0));
//    assert(a.nbytes() % sizeof(Input_type) == 0);
//    auto end_input = begin_input + (a.nbytes() / sizeof(Input_type));
//    std::transform(begin_input, end_input, result.mutable_data(0), F);
//    return result;
//}
//
//auto vertices_buffer() {
//    //std::cerr << "Size: " << c3t3.triangulation().number_of_vertices() << " vs " << vertices.size() << std::endl;
//    ////std::cerr << "Size: " << c3t3.triangulation().number_of_vertices() << " vs " << vertices_tmp.size() << std::endl;
//    return Vertices_buffer{ ArrayWrapper::wrap(vertices.front().data(), vertices.size()) };
//    //return Vertices_buffer{ ArrayWrapper::wrap(vertices_tmp) };
//}
//auto cells_buffer() {
//    //std::cerr << "Cell size: " << c3t3.number_of_cells() << " vs " << cells.size() << std::endl;
//    return Cells_buffer{ ArrayWrapper::wrap(cells.front().data(), cells.size()) };
//}
//auto domain_tags_buffer() {
//    return Domain_tag_buffer{ ArrayWrapper::wrap(domain_tags.data(), domain_tags.size()) };
//}
//auto facets_buffer() {
//    return Facets_buffer{ ArrayWrapper::wrap(facets.front().data(), facets.size()) };
//}
//auto facet_tags_buffer() {
//    return Facet_tag_buffer{ ArrayWrapper::wrap(facet_tags.front().data(), facet_tags.size()) };
//}
//auto squared_distances(Id_type f1, Id_type f2, py::array_t<double, py::array::c_style> pts) {
//    auto distance_functor = std::function<double(const Pt&)>{ squared_distance_functor(f1, f2) };
//    return apply(distance_functor, pts);
//}
//

//int compute_distances(const std::string& filename, const C3t3& c3t3)
//{
//
//	std::cout << "Extracting point data from file: " << filename << std::endl;
//	std::ifstream file(filename, std::ios_base::in);
//	if (!file.good()) {
//		file.close();
//		std::cerr << "Something went wrong in extracting points." << std::endl;
//		return -1;
//	}
//	else {
//		const DTree dtree(c3t3);
//		typedef C3t3::Surface_patch_index SPI;
//		std::vector<SPI> patch_indexes;
//		dtree.patch_indexes(std::back_inserter(patch_indexes));
//		std::ofstream index_dictionnary("index_dictionnary.txt");
//		for (std::vector<SPI>::const_iterator p = patch_indexes.begin(); p != patch_indexes.end(); ++p) {
//			index_dictionnary << "Index: " << p->first << "," << p->second << std::endl;
//		}
//		index_dictionnary.close();
//		std::ofstream patch_distances("distances_to_patches.txt");
//		DTree::Point P;
//		while (file.good()) {
//			file >> P;
//			if (!file.good()) break;
//			patch_distances << P;
//			for (std::vector<SPI>::const_iterator p = patch_indexes.begin(); p != patch_indexes.end(); ++p) {
//				patch_distances << " " << dtree.squared_distance(P, *p);
//			}
//			patch_distances << std::endl;
//		}
//	}
//	file.close();
//	return 0;
//
//}

//C++17 !
// #include <functional>
//int func(double) { return 0; }
//int main() {
//	std::function f{ func }; // guide #1 deduces function<int(double)>
//	int i = 5;
//	std::function g = [&](double) { return i; }; // guide #2 deduces function<int(double)>
//}

