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
#include <cstddef>
#include <iostream>
struct PartitioningInformationWrapper {
   typedef int integral_type;
   const integral_type* rowl_to_rowg;
   const integral_type* coll_to_colg;
   const integral_type nb_well_inj_local;
   const integral_type nb_well_prod_local;
   const integral_type nb_node_own;
   const integral_type nb_frac_own;
   const integral_type nb_node_local;
   const integral_type nb_frac_local;
   const integral_type nb_comp_thermique;
   const integral_type nb_well_inj_own;
   const integral_type nb_well_prod_own;
   PartitioningInformationWrapper()
       : rowl_to_rowg{nullptr},
         coll_to_colg{nullptr},
         nb_well_inj_local{0},
         nb_well_prod_local{0},
         nb_node_own{0},
         nb_frac_own{0},
         nb_node_local{0},
         nb_frac_local{0},
         nb_comp_thermique{0},
         nb_well_inj_own{0},
         nb_well_prod_own{0} {}
   auto nb_rowl() const {
      assert(rowl_to_rowg);
      return static_cast<std::size_t>(nb_node_own + nb_frac_own +
                                      nb_well_inj_own + nb_well_prod_own);
   }
   auto nb_coll() const {
      assert(coll_to_colg);
      return static_cast<std::size_t>(nb_node_local + nb_frac_local +
                                      nb_well_inj_local + nb_well_prod_local);
   }
};
