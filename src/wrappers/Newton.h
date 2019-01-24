//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include "NewtonIncrements.h"

struct LinearSolver
{
    double tolerance;
    std::size_t maximum_number_of_iterations;
    std::size_t failures;
    std::size_t total_linear_solver_iterations;
}

struct Newton
{
// RW
    double tolerance;
    std::size_t maximum_number_of_iterations;
//RO
    std::vector<double> residuals;
    std::size_t total_iterations;
    std::size_t failures;
    NewtonIncrements increments;
    void reset(tol, std::size_t int)
};
