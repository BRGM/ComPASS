#pragma once

#include <array>

template <std::size_t nb_components, std::size_t nb_phases,
          typename Real_type = double>
struct Model {
   using Context = int;
   using Real = Real_type;
   static constexpr std::size_t nc = nb_components;
   static constexpr std::size_t np = nb_phases;
   // FIXME: preprocessor directives to be removed!
#ifdef _THERMIQUE_
   static constexpr std::size_t nbdof = nc + 1;
#else
   static constexpr std::size_t nbdof = nc;
#endif
   using Component_vector = std::array<Real, nc>;
   using Phase_component_matrix = std::array<Component_vector, np>;
   using Phase_vector = std::array<Real, np>;
   using Accumulation_vector = std::array<Real, nbdof>;
};
