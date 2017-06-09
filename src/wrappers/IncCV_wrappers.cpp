#include <cstddef>
#include <array>

#include "XArrayWrapper.h"

template <std::size_t nb_components, std::size_t nb_phases>
struct Model {
	static constexpr std::size_t nc = nb_components;
	static constexpr std::size_t np = nb_phases;
};

template <typename Model_type, typename Real_type = double>
struct IncCV {
	typedef int Context;
	typedef Real_type Real;
	static constexpr std::size_t nc = Model_type::nc;
	static constexpr std::size_t np = Model_type::np;
	// FIXME: preprocessor directives to be removed!
#ifdef _THERMIQUE_
	static constexpr std::size_t nbdof = nc + 1;
#else
	static constexpr std::size_t nbdof = nc;
#endif
	typedef std::array<Real, nc> Component_vector;
	typedef std::array<Real, np> Phase_vector;
	typedef std::array<Component_vector, np> Phase_component_matrix;
	typedef std::array<Real, nbdof> Accumulation_vector;
	Context context;
	Real p;
	Real T;
	Phase_component_matrix C;
	Phase_vector S;
	Accumulation_vector accumulation;
};

// FIXME: This is to be removed later
typedef Model<1, 2> Model_1c2p;
typedef IncCV<Model_1c2p> X1c2p;
typedef X1c2p X;

typedef XArrayWrapper<X> StateArray;

// Fortran functions
extern "C"
{
	int number_of_components();
	int number_of_phases();
	int size_of_unknowns();
	void retrieve_dirichlet_node_states(StateArray&);
	void retrieve_node_states(StateArray&);
	void retrieve_fracture_states(StateArray&);
	void retrieve_cell_states(StateArray&);
}

#include "IncCV_wrappers.h"
#include <pybind11/numpy.h>
//#include <pybind11/stl.h>

template<typename AttributeType, typename PyClass>
auto add_attribute_array(PyClass& states, const char *name, std::size_t offset,
	std::vector<std::size_t> shape = std::vector<std::size_t>{}, std::vector<std::size_t> strides = std::vector<std::size_t>{}) {
	states.def_property_readonly(name, [=](py::object& self) {
		auto wrapper = self.cast<StateArray*>(); // C++ pointer to underlying StateArray instance
		auto attribute_position = reinterpret_cast<const AttributeType*>(reinterpret_cast<const unsigned char*>(wrapper->pointer) + offset);
		auto final_shape = std::vector<std::size_t>{ { wrapper->length } };
		std::copy(shape.begin(), shape.end(), std::back_inserter(final_shape));
		auto final_strides = std::vector<std::size_t>{ { sizeof(X) } };
		std::copy(strides.begin(), strides.end(), std::back_inserter(final_strides));
		return py::array_t<AttributeType, py::array::c_style>{ final_shape, final_strides, attribute_position, self };
	},
		py::return_value_policy::reference_internal // py::keep_alive<0, 1>(): because the StateArray instance must be kept alive as long as the attribute is used
		);
};

void add_IncCV_wrappers(py::module& module)
{

	module.def("check_IncCV",
		[]() {
		py::print("Check IncCV:", number_of_components(), number_of_phases(), size_of_unknowns(), sizeof(X));
	});

	// FUTURE: all the following until next FUTURE tag shall be useless soon (cf infra)
	auto PyStateArray = py::class_<StateArray>(module, "States");
	add_attribute_array<X::Context>(PyStateArray, "context", offsetof(X, context));
	add_attribute_array<X::Real>(PyStateArray, "p", offsetof(X, p));
	add_attribute_array<X::Real>(PyStateArray, "T", offsetof(X, T));
	add_attribute_array<X::Real>(PyStateArray, "C", offsetof(X, C), { X::np, X::nc }, { X::nc * sizeof(X::Real), sizeof(X::Real) });
	add_attribute_array<X::Real>(PyStateArray, "S", offsetof(X, S), { X::np }, { sizeof(X::Real) });
	add_attribute_array<X::Real>(PyStateArray, "accumulation", offsetof(X, accumulation), { X::nbdof }, { sizeof(X::Real) });

	module.def("dirichlet_node_states", []() { return StateArray::retrieve(retrieve_dirichlet_node_states); });
	module.def("node_states", []() { return StateArray::retrieve(retrieve_node_states); });
	module.def("fracture_states", []() { return StateArray::retrieve(retrieve_fracture_states); });
	module.def("cell_states", []() { return StateArray::retrieve(retrieve_cell_states); });

	// FUTURE: PYBIND11_NUMPY_DTYPE shall support arrays in forthcoming release of pybind11
	// PYBIND11_NUMPY_DTYPE(X, context, p, T, C, S, accumulation);
	// Yet with a slightly cumbersome syntax: X['p'] to access pressure instead of X.p
	//PYBIND11_NUMPY_DTYPE(X, context, p, T);
	//module.def("node_states", []() {
	//	auto wrapper = StateArray{};
	//	retrieve_boundary_states(wrapper);
	//	return py::array_t<X, py::array::c_style>{wrapper.length, wrapper.pointer};
	//});

}
