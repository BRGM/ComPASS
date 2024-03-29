//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <cassert>
#include <vector>

extern "C" {
std::size_t nb_primary_variables();
std::size_t nb_nodes();
std::size_t nb_fractures();
std::size_t nb_cells();
std::size_t nb_injectors();
std::size_t nb_producers();
std::size_t nb_mswell_nodes();
}

struct NewtonIncrements {
   template <typename T>
   struct Pointers {
      T* nodes;
      T* fractures;
      T* cells;
      T* injectors;
      T* producers;
      T* mswell_nodes;
   };
   std::vector<double> nodes;
   std::vector<double> fractures;
   std::vector<double> cells;
   std::vector<double> injectors;
   std::vector<double> producers;
   std::vector<double> mswell_nodes;

   NewtonIncrements() = default;
   NewtonIncrements(const NewtonIncrements&) = default;
   NewtonIncrements(NewtonIncrements&&) = default;
   static constexpr std::size_t npv() {
// FIXME: should take into account the presence of components in phases (boolean
// table)
// FIXME: should coincide with DefModel
#ifdef _THERMIQUE_
#ifndef _WITH_FREEFLOW_STRUCTURES_
      constexpr std::size_t n =
          2 + (1 + ComPASS_NUMBER_OF_COMPONENTS) * ComPASS_NUMBER_OF_PHASES;
#else
      constexpr std::size_t n =
          2 + (2 + ComPASS_NUMBER_OF_COMPONENTS) * ComPASS_NUMBER_OF_PHASES;
#endif
#else   // _THERMIQUE_
      constexpr std::size_t n =
          1 + (1 + ComPASS_NUMBER_OF_COMPONENTS) * ComPASS_NUMBER_OF_PHASES;
#endif  // _THERMIQUE_
      return n;
   }
   void init() {
      // CHECKME: fill operations should not be necessary
      nodes.resize(npv() * nb_nodes());
      std::fill(begin(nodes), end(nodes), 0);
      fractures.resize(npv() * nb_fractures());
      std::fill(begin(fractures), end(fractures), 0);
      cells.resize(npv() * nb_cells());
      std::fill(begin(cells), end(cells), 0);
      injectors.resize(nb_injectors());
      std::fill(begin(injectors), end(injectors), 0);
      producers.resize(nb_producers());
      std::fill(begin(producers), end(producers), 0);
      mswell_nodes.resize(npv() * nb_mswell_nodes());
      std::fill(begin(mswell_nodes), end(mswell_nodes), 0);
   }
   auto pointers() {
      return Pointers<double>{nodes.data(),     fractures.data(),
                              cells.data(),     injectors.data(),
                              producers.data(), mswell_nodes.data()};
   }
   auto pointers() const {
      return Pointers<const double>{nodes.data(),     fractures.data(),
                                    cells.data(),     injectors.data(),
                                    producers.data(), mswell_nodes.data()};
   }
};
