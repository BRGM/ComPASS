#pragma once

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
    typedef std::array<Component_vector, np> Phase_component_matrix;
    typedef std::array<Real, np> Phase_vector;
    typedef std::array<Real, nbdof> Accumulation_vector;
    Context context;
    Real p;
    Real T;
    Phase_component_matrix C;
    Phase_vector S;
    Accumulation_vector accumulation;
#ifdef _WIP_FREEFLOW_STRUCTURES_
    Phase_vector FreeFlow_phase_flowrate; // molar flowrate in the freeflow (atmosphere) at the interface
#endif // _WIP_FREEFLOW_STRUCTURES_
};

template <typename Model_type, typename Real_type = double>
struct NeumannBoundaryConditions {
    typedef Real_type Real;
    static constexpr std::size_t nc = Model_type::nc;
    typedef std::array<Real, nc> Component_vector;
    Component_vector molar_flux;
    Real heat_flux;
};

// FIXME: This is to be removed later
typedef IncCV<Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>> X;
typedef XArrayWrapper<X> StateArray;
using NeumannBC = NeumannBoundaryConditions<Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>>;
//typedef XArrayWrapper<Neumann> NeumannArray;
