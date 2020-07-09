#pragma once

template <typename Model_type, typename Real_type = double>
struct NeumannBoundaryConditions {
   using Model = Model_type;
   typename Model::Component_vector molar_flux;
   typename Model::Real heat_flux;
};

#include "Model.h"

using NeumannBC = NeumannBoundaryConditions<
    Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>>;
