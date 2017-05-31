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
	typedef std::array<Real, nc> Component_vector;
	typedef std::array<Real, np> Phase_vector;
	typedef std::array<Component_vector, np> Phase_component_matrix;
	// FIXME: preprocessor directives to be removed!
#ifdef _THERMIQUE_
	typedef std::array<Real, nc + 1> Accumulation_vector;
#else
	typedef std::array<Real, nc> Accumulation_vector;
#endif
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
	void retrieve_boundary_node_states(StateArray&);
	void retrieve_node_states(StateArray&);
	void retrieve_fracture_states(StateArray&);
	void retrieve_cell_states(StateArray&);
}

#include "IncCV_wrappers.h"
#include <pybind11/numpy.h>

void add_IncCV_wrappers(py::module& module)
{

	module.def("check_IncCV",
		[]() {
		py::print("Check IncCV:", number_of_components(), number_of_phases(), size_of_unknowns(), sizeof(X));
	});

	// FUTURE: all the following until next FUTURE tag shall be useless soon (cf infra)
	py::class_<StateArray>(module, "States")
		.def_property_readonly("context", [](const StateArray& wrapper) {
		return py::array_t<X::Context, py::array::c_style>{
			{wrapper.length}, { sizeof(X) }, &(wrapper.pointer->context)
		};
	})
		.def_property_readonly("p", [](const StateArray& wrapper) {
		return py::array_t<X::Real, py::array::c_style>{
			{wrapper.length}, { sizeof(X) }, &(wrapper.pointer->p)
		};
	})
		.def_property_readonly("T", [](const StateArray& wrapper) {
		return py::array_t<X::Real, py::array::c_style>{
			{wrapper.length}, { sizeof(X) }, &(wrapper.pointer->T)
		};
	})
		.def_property_readonly("C", [](const StateArray& wrapper) {
		return py::array_t<X::Real, py::array::c_style>{
			{wrapper.length, X::np, X::nc},
			{ wrapper.length * sizeof(X), X::nc * sizeof(X::Real), sizeof(X::Real) },
				wrapper.pointer->C.front().data()
		};
	})
		.def_property_readonly("S", [](const StateArray& wrapper) {
		return py::array_t<X::Real, py::array::c_style>{
			{wrapper.length, X::np},
			{ wrapper.length * sizeof(X), sizeof(X::Real) },
				wrapper.pointer->S.data()
		};
	})
		.def_property_readonly("accumulation", [](const StateArray& wrapper) {
		return py::array_t<X::Real, py::array::c_style>{
			// CHECKME: We use std:array::size member here because it depends on _THERMIQUE_ directive
			{wrapper.length, wrapper.pointer->accumulation.size()},
			{ wrapper.length * sizeof(X), sizeof(X::Real) },
				wrapper.pointer->accumulation.data()
		};
	});

	module.def("boundary_node_states", []() { return StateArray::retrieve(retrieve_boundary_node_states); });
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
