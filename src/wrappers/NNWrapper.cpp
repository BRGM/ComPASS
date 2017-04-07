#include <cstddef>
#include <cassert>
#include <iterator>

#include<iostream>

#include <StringWrapper.h>

#include <pybind11/pybind11.h>

struct Vertices
{
	double * p;
	std::size_t nb_points;
	Vertices() :
		p(NULL),
		nb_points(0)
	{}
};

struct COC_container
{
protected:
	int * begin_;
	int * end_;
	COC_container() = delete;
	COC_container(int * content, int length) :
		begin_{ content },
		end_{ content + length }
	{
		assert(length >= 0);
	}
public:
	int * begin() const { return begin_; }
	int * end() const { return end_; }
	std::size_t length() const { return std::distance(begin_, end_); }
	friend class COC_iterator;
};

class COC_iterator
{
protected:
	// This is const as COC structure is not supposed to be changed
	const int * container_offset;
	int * container_content;
	COC_iterator() = delete;
	COC_iterator(int * offset, int * content) :
		container_offset{ offset },
		container_content{ content }
	{}
	auto container_length() const {
		return *std::next(container_offset) - *container_offset;
	}
public:
	auto operator*() const {
		return COC_container{ container_content, container_length() };
	}
	auto operator<(const COC_iterator& other) const {
		return container_content < other.container_content;
	}
	auto operator==(const COC_iterator& other) const {
		return container_content == other.container_content;
	}
	auto operator!=(const COC_iterator& other) const {
		return container_content != other.container_content;
	}
	COC_iterator& operator++() {
		// WARNING: container_content must be advanced before offset
		//          otherwise length is unvalidated
		std::advance(container_content, container_length());
		std::advance(container_offset, 1);
	}
	friend class COC;
};

/** Container of containers. */
class COC
{
protected:
	int nb_containers;
	int * container_offset;
	int * container_content;
public:
	auto begin() const {
		return COC_iterator{ container_offset, container_content };
	}
	auto end() const {
		assert(nb_containers >= 0);
		assert(container_offset[nb_containers] >= 0);
		auto offset = container_offset;
		auto content = container_content;
		std::advance(offset, nb_containers);
		std::advance(content, container_offset[nb_containers]);
		return COC_iterator{ offset, content };
	}
	auto operator[](const int i) {
		assert(i >= 0);
		assert(i < nb_containers);
		return *COC_iterator{ &container_offset[i], &container_content[container_offset[i]] };
	}
};

struct MeshConnectivity
{
	COC NodebyCell;
	COC NodebyFace;
	COC FacebyCell;
	COC CellbyNode;
	COC CellbyFace;
	COC CellbyCell;
};

std::ostream& operator<<(std::ostream& os, const MeshConnectivity& connectivity)
{

	// Classical nested loops
	os << "Mesh cell nodes:" << std::endl;
	for (auto && cell : connectivity.NodebyCell) {
		for (auto && node : cell) {
			os << node << " ";
		}
		os << std::endl;
	}
	// Alternative way of dumping nodes using STL facilities
	os << "Mesh cell faces:" << std::endl;
	for (auto && face : connectivity.NodebyFace) {
		std::copy(std::begin(face), std::end(face), std::ostream_iterator<int>(os, " "));
		os << std::endl;
	}
	return os;

}

extern "C"
{
	void NN_init(const StringWrapper&, const StringWrapper&, const StringWrapper&);
	void NN_init_warmup_and_read_mesh(const StringWrapper&, const StringWrapper&);
	void NN_init_warmup(const StringWrapper&);
	void NN_init_read_mesh(const StringWrapper&);
	void NN_init_build_grid(double, double, double, double, double, double, int, int, int);
	void NN_init_phase2(const StringWrapper&);
	void NN_main(int, const StringWrapper&);
	void NN_finalize();
	void retrieve_vertices(Vertices&);
	void retrieve_mesh_connectivity(MeshConnectivity&);
}

namespace py = pybind11;

PYBIND11_PLUGIN(ComPASS)
{

	py::module module("ComPASS", "pybind11 ComPASS library interface");

	module.def("init",
		[](const std::string& MeshFile, const std::string& LogFile, const std::string& OutputDir) {
		NN_init(MeshFile, LogFile, OutputDir);
	},
		"Initialisation of ComPASS.");

	module.def("init_warmup_and_read_mesh",
	[](const std::string& MeshFile, const std::string& LogFile) {
		NN_init_warmup_and_read_mesh(MeshFile, LogFile);
	},
	"Initialisation of ComPASS up to the point where the global mesh is loaded.  Next logical step is init_phase2.");

	module.def("init_warmup",
	[](const std::string& LogFile) {
		NN_init_warmup(LogFile);
	},
	"Initialisation of ComPASS - warmup phase. Next logical step is init_read_mesh or init_build_grid.");

	module.def("init_read_mesh",
	[](const std::string& MeshFile) {
		NN_init_read_mesh(MeshFile);
	},
	"Initialisation of ComPASS - read the mesh file. Next logical step is init_phase2.");

	//module.def("init_build_grid", [](py::sequence origin, py::sequence extent, py::sequence shape) {
	//	NN_init_build_grid(
	//		origin[0], origin[1], origin[2], 
	//		extent[0], extent[1], extent[2], 
	//		shape[0], shape[1], shape[2] 
	//	);
	//},
	//"Initialisation of ComPASS - build a cartesian grid. Next logical step is init_phase2.");
	//module.def("init_build_grid", [](double Ox, double Oy, double Oz, double lx, double ly, double lz, int nx, int ny, int nz) {
	//	NN_init_build_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz);
	//},
			   //"Initialisation of ComPASS - build a cartesian grid. Next logical step is init_phase2.");
	module.def("init_build_grid", &NN_init_build_grid,
			   "Initialisation of ComPASS - build a cartesian grid. Next logical step is init_phase2.");
	
	module.def("init_phase2",
		[](const std::string& OutputDir) {
		NN_init_phase2(OutputDir);
	},
		"Initialisation of ComPASS phase 2 : partition and distribute.");
	
	module.def("main_loop", [](int TimeIter, const std::string& OutputDir) { NN_main(TimeIter, OutputDir); },
		"Main loop of ComPASS.");
	
	module.def("finalize", &NN_finalize, "Cleans ComPASS data structures.");

	py::class_<Vertices>(module, "VerticesHandle", py::buffer_protocol())
		.def_buffer([](Vertices &V) -> py::buffer_info {
		return py::buffer_info(
			V.p,                               /* Pointer to buffer */
			sizeof(double),                          /* Size of one scalar */
			py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
			2,                                      /* Number of dimensions */
			{ V.nb_points, 3 },                 /* Buffer dimensions */
			{ sizeof(double) * 3,             /* Strides (in bytes) for each index */
			  sizeof(double) }
		);
	});
	
	module.def("get_vertices", []() {
		Vertices V;
		retrieve_vertices(V);
		return V;
	},
		"Get node coordinates."
		);

	py::class_<COC_container>(module, "COCcontainer", py::buffer_protocol())
		.def_buffer([](COC_container &cocc) -> py::buffer_info {
		return py::buffer_info(
			std::begin(cocc),                               /* Pointer to buffer */
			sizeof(int),                          /* Size of one scalar */
			py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
			1,                                      /* Number of dimensions */
			{ cocc.length() },                 /* Buffer dimensions */
			{ sizeof(int) }
		);
	})
		.def("__len__", [](const COC_container &cocc) { return cocc.length(); });

	py::class_<COC_iterator>(module, "COCiterator");

	py::class_<COC>(module, "COC")
		.def("__iter__", [](COC& coc) {
		return py::make_iterator(coc.begin(), coc.end());
	},
			py::keep_alive<0, 1>() /* Keep COC alive while iterator is used */
		)
	.def("__getitem__", [](COC& coc, int i) -> COC_container { return coc[i]; })
		;

	py::class_<MeshConnectivity>(module, "MeshConnectivity")
		.def_readwrite("NodebyCell", &MeshConnectivity::NodebyCell)
		.def_readwrite("NodebyFace", &MeshConnectivity::NodebyFace)
		.def_readwrite("FacebyCell", &MeshConnectivity::FacebyCell)
		.def_readwrite("CellbyNode", &MeshConnectivity::CellbyNode)
		.def_readwrite("CellbyFace", &MeshConnectivity::CellbyFace)
		.def_readwrite("CellbyCell", &MeshConnectivity::CellbyCell);

	module.def("get_connectivity", []() {
		MeshConnectivity connectivity;
		retrieve_mesh_connectivity(connectivity);
		//std::cout << connectivity << std::endl;
		return connectivity;
	},
		"Get mesh connectivity."
		);

	return module.ptr();

}
