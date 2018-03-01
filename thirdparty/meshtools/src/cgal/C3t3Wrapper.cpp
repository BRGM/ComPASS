#include <iostream>
#include <fstream>
#include <utility>
#include <iterator>
#include <memory>
#include <functional>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Mesh_domain_with_polyline_features_3.h"
#include "CGAL/Mesh_triangulation_3.h"
#include "CGAL/Mesh_complex_3_in_triangulation_3.h"

#include "DummyDomain.h"
#include "DistanceTrees.h"
#include "ArrayWrapper.h"
#include "PyBuffer_wrappers.h"


typedef CGAL::Epick Kernel;
typedef int SubdomainIndex;
typedef std::pair<SubdomainIndex, SubdomainIndex> SurfacePatchIndex;

typedef ExperimentalOracle::DummyDomain<Kernel, SubdomainIndex, SurfacePatchIndex> Domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Domain> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain, Mesh_domain::R, CGAL::Sequential_tag>::type Triangulation;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Triangulation> C3t3;
typedef ExperimentalBRGM::DistanceTrees<C3t3> DTree;
typedef DTree::Point DTP;

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

//C++17 !
// #include <functional>
//int func(double) { return 0; }
//int main() {
//	std::function f{ func }; // guide #1 deduces function<int(double)>
//	int i = 5;
//	std::function g = [&](double) { return i; }; // guide #2 deduces function<int(double)>
//}

template<typename Result, typename Input>
auto apply(std::function<Result(Input)> F, py::array_t<double, py::array::c_style> a)
{
	typedef typename std::remove_reference<Input>::type Input_type;
	typedef py::array_t<Result, py::array::c_style> Result_array;
	if (a.ndim() == 0) return Result_array{ {} };
	std::vector<std::size_t> result_shape;
	// IMPROVE: set as boolean template
	if (a.itemsize() == sizeof(Input_type)) {
        auto shape = a.request().shape;
        result_shape.reserve(shape.size());
        for (auto&& n : shape) {
            result_shape.push_back(n);
        }
	}
	else if (a.shape(a.ndim() - 1) * a.itemsize() == sizeof(Input_type)) {
		const auto info = a.request();
		result_shape = std::vector<std::size_t>{ info.shape.begin(), info.shape.begin() + a.ndim() - 1 };
	}
	else {
		throw std::runtime_error("Incompatible size in memory!");
	}
	auto result = Result_array{ result_shape };
	// WARNING: If Input_type is not const cast will fail as a.data is const: use a.mutable_data
	auto begin_input = reinterpret_cast<Input_type *>(a.data(0));
	assert(a.nbytes() % sizeof(Input_type) == 0);
	auto end_input = begin_input + (a.nbytes() / sizeof(Input_type));
	std::transform(begin_input, end_input, result.mutable_data(0), F);
	return result;
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
	typedef PyBufferWrapper<double, dim> Vertices_buffer;
	typedef PyBufferWrapper<Id_type, dim + 1> Cells_buffer;
	typedef PyBufferWrapper<Id_type, dim> Facets_buffer;
	typedef PyBufferWrapper<Id_type> Domain_tag_buffer;
	typedef PyBufferWrapper<Id_type, 2> Facet_tag_buffer;
protected:
	C3t3 c3t3;
	//std::vector<Field_type> vertices_tmp;
	std::vector<Vertex> vertices;
	std::vector<Simplex> cells;
	std::vector<Facet> facets;
	std::vector<Domain_tag> domain_tags;
	std::vector<Facet_tag> facet_tags;
	std::shared_ptr<DTree> dtree;
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
		dtree = std::make_shared<DTree>(complex);
	}
public:
	void reload(const std::string& filename) {
		load_c3t3(filename, c3t3);
		update_info();
	}
	auto vertices_buffer() {
		//std::cerr << "Size: " << c3t3.triangulation().number_of_vertices() << " vs " << vertices.size() << std::endl;
		////std::cerr << "Size: " << c3t3.triangulation().number_of_vertices() << " vs " << vertices_tmp.size() << std::endl;
		return Vertices_buffer{ ArrayWrapper::wrap(vertices.front().data(), vertices.size()) };
		//return Vertices_buffer{ ArrayWrapper::wrap(vertices_tmp) };
	}
	auto cells_buffer() {
		//std::cerr << "Cell size: " << c3t3.number_of_cells() << " vs " << cells.size() << std::endl;
		return Cells_buffer{ ArrayWrapper::wrap(cells.front().data(), cells.size()) };
	}
	auto domain_tags_buffer() {
		return Domain_tag_buffer{ ArrayWrapper::wrap(domain_tags.data(), domain_tags.size()) };
	}
	auto facets_buffer() {
		return Facets_buffer{ ArrayWrapper::wrap(facets.front().data(), facets.size()) };
	}
	auto facet_tags_buffer() {
		return Facet_tag_buffer{ ArrayWrapper::wrap(facet_tags.front().data(), facet_tags.size()) };
	}
	auto squared_distance_functor(Id_type f1, Id_type f2) const {
		auto tree = std::shared_ptr<DTree>{ dtree };
		auto spi = SurfacePatchIndex{ f1, f2 };
		return [tree, spi](const Pt& P) -> double {
			static_assert(sizeof(DTree::Point) == sizeof(Pt), "Incompatible Point types !");
			return tree->squared_distance(reinterpret_cast<const DTree::Point&>(P), spi);
		};
	}
	auto squared_distances(Id_type f1, Id_type f2, py::array_t<double, py::array::c_style> pts) {
		auto distance_functor = std::function<double(const Pt&)>{ squared_distance_functor(f1, f2) };
		return apply(distance_functor, pts);
	}
};

struct Foo
{
	py::array_t<double, py::array::c_style> a;
	py::array_t<double, py::array::c_style> b;
	py::array_t<int, py::array::c_style> c;
	Foo() :
		a{ { 2, 3 } },
		b{ { 5, 2 } },
		c{ { 4, 4 } } {
		// Hack to make the array const from: https://github.com/pybind/pybind11/issues/481
		reinterpret_cast<py::detail::PyArray_Proxy*>(c.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
	}
};

void f(py::array_t<double, py::array::c_style> a)
{
	auto info = a.request();
	const auto ndim = info.ndim;
	if (ndim < 1 || ndim > 2)
		throw std::runtime_error("Wrong number of dimensions!");
	assert(shape.size() == ndim);
	if(info.shape.back() != 3)
		throw std::runtime_error("3D points only!");
	const std::size_t n = info.shape[0];
	std::size_t nbpts = ndim == 1 ? 1 : info.shape[0];
	double *p = static_cast<double *>(info.ptr);
	for(;nbpts!=0; --nbpts) {
		std::cerr << *p << " ";
		++p;
		std::cerr << *p << " ";
		++p;
		std::cerr << *p << std::endl;
		++p;
	}
}

void g(py::array_t<double, py::array::c_style> a)
{

	std::cerr << "g ndim:" << a.ndim() << std::endl;
	std::cerr << "g shape:" << a.shape(0) << " " << a.shape(1) << std::endl;
	std::cerr << a.data(0, 0) << " " << a.data(0, 1) << std::endl;
	//const double *p = a.data(1);
	//const Pt& P = *(static_cast<const Pt*>(reinterpret_cast<const void *>(p)));
	const Pt& P = *(reinterpret_cast<const Pt *>(a.data(1)));
	std::cerr << "P:" << P.x << " " << P.y << std::endl;
}

double f1(const Pt& P)
{
	return P.x + P.y + P.z;
}

int f2(const Pt& P)
{
	return static_cast<int>(P.x + P.y + P.z);
}

PYBIND11_MODULE(C3t3Wrapper, module)
{

    module.doc() = "pybind11 homemade CGAL C3t3 interface";

	C3t3Wrapper::Vertices_buffer::add_buffer_class(module, "Vertices");
	C3t3Wrapper::Cells_buffer::add_buffer_class(module, "Cells");
	C3t3Wrapper::Domain_tag_buffer::add_buffer_class(module, "DomainTags");
	C3t3Wrapper::Facets_buffer::add_buffer_class(module, "Facets");
	C3t3Wrapper::Facet_tag_buffer::add_buffer_class(module, "FacetTags");

	py::class_<C3t3Wrapper>(module, "C3t3")
		.def(py::init<>())
		.def("reload", &C3t3Wrapper::reload)
		.def("vertices", &C3t3Wrapper::vertices_buffer, py::return_value_policy::reference_internal)
		.def("cells", &C3t3Wrapper::cells_buffer, py::return_value_policy::reference_internal)
		.def("domain_tags", &C3t3Wrapper::domain_tags_buffer, py::return_value_policy::reference_internal)
		.def("facets", &C3t3Wrapper::facets_buffer, py::return_value_policy::reference_internal)
		.def("facet_tags", &C3t3Wrapper::facet_tags_buffer, py::return_value_policy::reference_internal)
		.def("squared_distances", &C3t3Wrapper::squared_distances);

	//module.def("toto", f);

	//module.def("tutu", [](std::size_t n) {
	//	auto result = py::array_t<double, 2>({ n, 3 });
	//	return result;
	//});

	PYBIND11_NUMPY_DTYPE(Pt, x, y, z);

	//module.def("test", [](std::size_t n) {
	//	auto result = py::array_t<Pt, 2>({ n, 1 });
	//	return result;
	//});

	//module.def("f1", [](py::array_t<double, py::array::c_style> a) {
	//	// IMPROVE: C++17: type deduction on std::function, recast_function is no longer needed
	//	return apply(recast_function(f1), a);
	//});

	//module.def("f2", [](py::array_t<double, py::array::c_style> a) {
	//	// IMPROVE: C++17: type deduction on std::function, recast_function is no longer needed
	//	return apply(recast_function(f2), a);
	//});

	py::class_<Foo>(module, "Foo")
		.def(py::init<>())
		.def_readwrite("a", &Foo::a)
		.def_readonly("b", &Foo::b) // you can change the content of b but not the property : shape, dtype, allocated chunk of memory
		.def_readonly("c", &Foo::c); // you cannot change the content of c because it is const int

}

//int main(int argc, const char *argv[])
//{
//	C3t3Wrapper wrapper;
//	if (argc > 1) {
//		wrapper.reload(argv[1]);
//	}
//
//	//C3t3 c3t3;
//	//if(argc>1) {
//	//	load_c3t3(argv[1], c3t3);
//	//	if(argc>2) {
//	//		compute_distances(argv[2], c3t3);
//	//	}
//	//}
//	//else {
//	//	assert(argc > 0);
//	//	std::cerr << argv[0] << " c3t3file distances_results" << std::endl;
//	//}
//	return 0;
//
//}
