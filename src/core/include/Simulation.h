//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

namespace ComPASS {

// A dummy object used as singleton (one per proc)
struct Simulation {
   int proc_rank;
   Simulation(int p) : proc_rank(p) {}
};

void create_simulation_object(const int);
Simulation* simulation_object();

}  // namespace ComPASS
