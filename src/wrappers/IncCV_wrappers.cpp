//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <cstddef>
#include <array>
#include <sstream>

#include "StateObjects.h"

// Fortran functions
extern "C"
{
	int number_of_components();
	int number_of_phases();
	int size_of_unknowns();
	void retrieve_dirichlet_node_states(StateArray&);
	void retrieve_all_states(StateArray&);
	void retrieve_node_states(StateArray&);
	void retrieve_fracture_states(StateArray&);
	void retrieve_cell_states(StateArray&);
	void retrieve_own_dirichlet_node_states(StateArray&);
	void retrieve_own_node_states(StateArray&);
	void retrieve_own_fracture_states(StateArray&);
	void retrieve_own_cell_states(StateArray&);
	void retrieve_injection_whp(XArrayWrapper<double>&);
    void retrieve_production_whp(XArrayWrapper<double>&);
    void retrieve_production_whp(XArrayWrapper<double>&);
    void set_face_with_neumann_contribution(int, const NeumannBC&);
    void set_faces_with_neumann_contribution(int, const int*, const NeumannBC&);
    void set_face_neumann_contributions(int, int*, NeumannBC*);
    void set_fracture_edge_with_neumann_contribution(const int*, const NeumannBC&);
    void set_fracture_edges_with_neumann_contribution(int, const int*, const NeumannBC&);
    void set_fracture_edge_neumann_contributions(int, int*, NeumannBC*);
    void clear_all_neumann_contributions();
#ifndef NDEBUG
    void dump_incv_info();
#endif
}

#include "IncCV_wrappers.h"
#include <pybind11/numpy.h>
#include "PyXArrayWrapper.h"

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
		py::return_value_policy::reference_internal // py::keep_alive<0, 1>(): because the StateArray instance must be kept alive as long as the attribute is used - should be ok this is the default for def_property
		);
};

void add_IncCV_wrappers(py::module& module)
{

	module.def("check_IncCV",
		[]() {
		py::print("Check IncCV:", number_of_components(), number_of_phases(), size_of_unknowns(), sizeof(X));
	});

	py::class_<X>(module, "State")
        .def_readwrite("p", &X::p)
        .def_readwrite("T", &X::T)
        .def_property_readonly("S", [](py::object& self) {
            return py::array_t<X::Real, py::array::c_style>{ X::np, self.cast<X&>().S.data(), self};
        })
        .def_property_readonly("C", [](py::object& self) {
            return py::array_t<X::Real, py::array::c_style>{ {X::np, X::nc}, self.cast<X&>().C.data()->data(), self};
        })
        .def_property_readonly("accumulation",  [](py::object& self) {
            return py::array_t<X::Real, py::array::c_style>{X::nc, self.cast<X&>().accumulation.data(), self};
        })
        .def("__str__", [](const X& self) {
            std::basic_stringstream<char> s;
            s << "Context: " << static_cast<int>(self.context);
            s << " p=" << self.p;
            s << " T=" << self.T;
            s << " S=(";
            for(auto&& Sk: self.S) s << Sk <<", ";            
            s << ") C=(";
            for(auto&& Ck: self.C) {
                s << "(" ;
                for(auto&& Cki: Ck) s << Cki <<", ";   
                s << "), ";         
            }
            s << ")";
            return py::str{ s.str() };
        })
    ;

    // FUTURE: all the following until next FUTURE tag shall be useless soon (cf infra)
	auto PyStateArray = py::class_<StateArray>(module, "States")
		.def("size", [](const StateArray& self) { return self.length; })
		.def_property_readonly("shape", [](const StateArray& self) { 
            return py::make_tuple(self.length);
        })
        .def("__getitem__", &StateArray::operator[], py::return_value_policy::reference)
        .def("__iter__", [](StateArray& self) {
		    return py::make_iterator(
                self.pointer, self.pointer + self.length
            );
        })
        .def("fill", &StateArray::fill)
        .def("set", &StateArray::fill)
    ;

	add_attribute_array<typename X::Context>(PyStateArray, "context", offsetof(X, context));
	add_attribute_array<typename X::Real>(PyStateArray, "p", offsetof(X, p));
	add_attribute_array<typename X::Real>(PyStateArray, "T", offsetof(X, T));
	add_attribute_array<typename X::Real>(PyStateArray, "C", offsetof(X, C), { X::np, X::nc }, { X::nc * sizeof(X::Real), sizeof(X::Real) });
	add_attribute_array<typename X::Real>(PyStateArray, "S", offsetof(X, S), { X::np }, { sizeof(X::Real) });
	add_attribute_array<typename X::Real>(PyStateArray, "accumulation", offsetof(X, accumulation), { X::nbdof }, { sizeof(X::Real) });
#ifdef _WIP_FREEFLOW_STRUCTURES_
    add_attribute_array<typename X::Real>(PyStateArray, "FreeFlow_phase_flowrate", offsetof(X, FreeFlow_phase_flowrate), { X::np }, { sizeof(X::Real) });
#endif // _WIP_FREEFLOW_STRUCTURES_
    //auto PyNeumannArray = py::class_<NeumannArray>(module, "NeumannContributions")
    //    .def("size", [](const NeumannArray& self) { return self.length; })
    //    .def_property_readonly("shape", [](const NeumannArray& self) { return py::make_tuple(self.length); });
    //add_attribute_array<typename Neumann::molar_flux>(PyNeumannArray, "molar_flux", offsetof(Neumann, molar_flux), { Neumann::nc }, { sizeof(X::Real) });
    //add_attribute_array<typename Neumann::Real>(PyNeumannArray, "heat_flux", offsetof(Neumann, heat_flux));

	module.def("dirichlet_node_states", []() { return StateArray::retrieve(retrieve_dirichlet_node_states); });
	module.def("all_states", []() { return StateArray::retrieve(retrieve_all_states); },
    R"doc(
        Retrieve all nodes, fractures and cell states in a contiguous structure.
    )doc");
	module.def("node_states", []() { return StateArray::retrieve(retrieve_node_states); });
	module.def("fracture_states", []() { return StateArray::retrieve(retrieve_fracture_states); });
	module.def("cell_states", []() { return StateArray::retrieve(retrieve_cell_states); });
	module.def("own_dirichlet_node_states", []() { return StateArray::retrieve(retrieve_own_dirichlet_node_states); });
	module.def("own_node_states", []() { return StateArray::retrieve(retrieve_own_node_states); });
	module.def("own_fracture_states", []() { return StateArray::retrieve(retrieve_own_fracture_states); });
	module.def("own_cell_states", []() { return StateArray::retrieve(retrieve_own_cell_states); });

	module.def("number_of_components", []() { return ComPASS_NUMBER_OF_COMPONENTS; });
	module.def("number_of_phases", []() { return ComPASS_NUMBER_OF_PHASES; });

	// FUTURE: PYBIND11_NUMPY_DTYPE shall support arrays in forthcoming release of pybind11
	// PYBIND11_NUMPY_DTYPE(X, context, p, T, C, S, accumulation);
	// Yet with a slightly cumbersome syntax: X['p'] to access pressure instead of X.p
	//PYBIND11_NUMPY_DTYPE(X, context, p, T);
	//module.def("node_states", []() {
	//	auto wrapper = StateArray{};
	//	retrieve_boundary_states(wrapper);
	//	return py::array_t<X, py::array::c_style>{wrapper.length, wrapper.pointer};
	//});
    
    py::class_<NeumannBC>(module, "NeumannBC")
        .def(py::init<>())
        .def_readwrite("heat_flux", &NeumannBC::heat_flux)
            .def_property_readonly("molar_flux", [](py::object& self) {
        auto data = self.cast<NeumannBC*>()->molar_flux.data(); // C++ pointer to underlying NeumannBC instance
        return py::array_t<NeumannBC::Real, py::array::c_style>({ NeumannBC::nc }, { sizeof(NeumannBC::Real) }, data, self);
    }) // return value policy defaults to py::return_value_policy::reference_internal which is ok
        ;

    module.def("set_Neumann_faces", [](
        py::array_t<int, py::array::c_style | py::array::forcecast> faces,
        const NeumannBC& condition
        ) {
        if (faces.ndim() != 1 || faces.size()==0) return;
        auto raw_faces = faces.unchecked<1>();
        set_faces_with_neumann_contribution(
            raw_faces.shape(0), raw_faces.data(0), condition
        );
    });

    module.def("set_Neumann_fracture_edges", [](
        py::array_t<int, py::array::c_style | py::array::forcecast> edges,
        const NeumannBC& condition
        ) {
        if (edges.ndim() != 2 || edges.size() == 0) return;
        auto raw_edges = edges.unchecked<2>();
        set_fracture_edges_with_neumann_contribution(
            raw_edges.shape(0), raw_edges.data(0, 0), condition
        );
    });

    module.def("clear_all_neumann_contributions", &clear_all_neumann_contributions);

    //PYBIND11_NUMPY_DTYPE(NeumannBC, molar_flux, heat_flux);
    //module.def("neumann_conditions", [](std::size_t n) {
    //    return py::array_t<NeumannBC, py::array::c_style>{ n };
    //});

    //module.def("node_states", []() {
    //	auto wrapper = StateArray{};
    //	retrieve_boundary_states(wrapper);
    //	return py::array_t<X, py::array::c_style>{wrapper.length, wrapper.pointer};
    //});


	add_array_wrapper(module, "injection_whp", retrieve_injection_whp);
	add_array_wrapper(module, "production_whp", retrieve_production_whp);

#ifndef NDEBUG
    module.def("dump_incv_info", &dump_incv_info);
#endif

}
