#pragma once
#include <array>

template <typename Model_type>
struct Xa {
   using Model = Model_type;
   typename Model::Real p;
   typename Model::Real T;
   typename Model::Phase_component_matrix C;
};

template <typename Model_type>
struct IncCV {
   using Model = Model_type;
   typename Model::Context context;
   typename Model::Real p;
   typename Model::Phase_vector pa;
   typename Model::Real T;
   typename Model::Phase_component_matrix C;
   typename Model::Phase_vector S;
   typename Model::Accumulation_vector accumulation;
#ifdef _WITH_FREEFLOW_STRUCTURES_
   typename Model::Phase_vector
       FreeFlow_phase_flowrate;  // molar flowrate in the freeflow
                                 // (atmosphere) at the interface
#endif                           // _WITH_FREEFLOW_STRUCTURES_
};

template <typename Model_type>
struct FFfarfield {
   using Model = Model_type;
   typename Model::Real p;
   typename Model::Phase_vector T;
   typename Model::Phase_component_matrix C;
   typename Model::Phase_vector imposed_flux;
   typename Model::Phase_vector Hm;
   typename Model::Real HT;
};

// FIXME: This is to be removed later

#include "Model.h"
#include "XArrayWrapper.h"

using X = IncCV<Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>>;
using StateArray = XArrayWrapper<X>;
using XFF =
    FFfarfield<Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>>;
using StateFFArray = XArrayWrapper<XFF>;
using Xalpha =
    Xa<Model<ComPASS_NUMBER_OF_COMPONENTS, ComPASS_NUMBER_OF_PHASES>>;
