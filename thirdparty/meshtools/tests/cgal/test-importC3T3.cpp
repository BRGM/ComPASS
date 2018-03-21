#include <iostream>
#include <fstream>
#include <utility>
#include <iterator>

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Mesh_domain_with_polyline_features_3.h"
#include "CGAL/Mesh_triangulation_3.h"
#include "CGAL/Mesh_complex_3_in_triangulation_3.h"

#include "DummyDomain.h"
#include "DistanceTrees.h"
#include "PyBuffer_wrappers.h"

typedef CGAL::Epick Kernel;
typedef int SubdomainIndex;
typedef std::pair<SubdomainIndex, SubdomainIndex> SurfacePatchIndex;

typedef ExperimentalOracle::DummyDomain<Kernel, SubdomainIndex, SurfacePatchIndex> Domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Domain> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain, Mesh_domain::R, CGAL::Sequential_tag>::type Triangulation;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Triangulation> C3t3;
typedef ExperimentalBRGM::DistanceTrees<C3t3> DTree;

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

int load_c3t3(const std::string& filename, C3t3& c3t3)
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
		return -1;
	}
	file.close();
	return 0;

}

int compute_distances(const std::string& filename, const C3t3& c3t3)
{

	std::cout << "Extracting point data from file: " << filename << std::endl;
	std::ifstream file(filename, std::ios_base::in);
	if (!file.good()) {
		file.close();
		std::cerr << "Something went wrong in extracting points." << std::endl;
		return -1;
	}
	else {
		const DTree dtree(c3t3);
		typedef C3t3::Surface_patch_index SPI;
		std::vector<SPI> patch_indexes;
		dtree.patch_indexes(std::back_inserter(patch_indexes));
		std::ofstream index_dictionnary("index_dictionnary.txt");
		for (std::vector<SPI>::const_iterator p = patch_indexes.begin(); p != patch_indexes.end(); ++p) {
			index_dictionnary << "Index: " << p->first << "," << p->second << std::endl;
		}
		index_dictionnary.close();
		std::ofstream patch_distances("distances_to_patches.txt");
		DTree::Point P;
		while (file.good()) {
			file >> P;
			if (!file.good()) break;
			patch_distances << P;
			for (std::vector<SPI>::const_iterator p = patch_indexes.begin(); p != patch_indexes.end(); ++p) {
				patch_distances << " " << dtree.squared_distance(P, *p);
			}
			patch_distances << std::endl;
		}
	}
	file.close();
	return 0;

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
	typedef std::array<Id_type, 2> Facet_tag;
	typedef PyBufferWrapper<double, dim> Vertices_buffer;
	typedef PyBufferWrapper<Id_type, dim + 1> Cells_buffer;
	typedef PyBufferWrapper<Id_type, dim> Facets_buffer;
	typedef PyBufferWrapper<Id_type, dim> Facet_tag_buffer;
protected:
	C3t3 c3t3;
	std::vector<Vertex> vertices;
	std::vector<Simplex> cells;
	std::vector<Facet> facets;
	std::vector<Facet_tag> facet_tags;
	void update_info() {
		const C3t3& complex = c3t3;
		const Triangulation& triangulation = complex.triangulation();
		vertices.clear();
		vertices.reserve(triangulation.number_of_vertices());
		assert(vertices.empty());
		cells.clear();
		cells.reserve(complex.number_of_cells_in_complex());
		assert(cells.empty());
		facets.clear();
		facets.reserve(complex.number_of_facets_in_complex());
		assert(facets.empty());
		facet_tags.clear();
		facet_tags.reserve(complex.number_of_facets_in_complex());
		assert(facet_tags.empty());
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
	}
public:
	void reload(const std::string& filename) {
		load_c3t3(filename, c3t3);
		update_info();
	}
};

int main(int argc, const char *argv[])
{
	C3t3Wrapper wrapper;
	if (argc > 1) {
		wrapper.reload(argv[1]);
	}

	//C3t3 c3t3;
	//if(argc>1) {
	//	load_c3t3(argv[1], c3t3);
	//	if(argc>2) {
	//		compute_distances(argv[2], c3t3);
	//	}
	//}
	//else {
	//	assert(argc > 0);
	//	std::cerr << argv[0] << " c3t3file distances_results" << std::endl;
	//}
	return 0;

}
