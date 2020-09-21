#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
from ..mpi import master_print
from ..options import get

IterativeSolverSettings = namedtuple(
    "IterativeSolverSettings",
    ["tolerance", "max_iterations", "restart_size"],
    defaults=[None, None, None],
)


class LinearSolver:
    """
    A base structure to hold the common parameters of ComPASS linear solvers
    """

    def __init__(self, linear_system):
        """
        :param linear_system: the linear system structure
        """
        self.failures = 0
        self.linear_system = linear_system
        if get("--linear_solver_view", False):
            master_print(self)

    def __str__(self):
        return "LinearSolver object view:"


class DirectSolver(LinearSolver):
    def __init__(self, linear_system):
        super().__init__(linear_system)

    def __str__(self):
        return f"{super().__str__()}\n   Direct"


class IterativeSolver(LinearSolver):
    def __init__(self, linear_system, settings):
        """
        :param settings: An IterativeSolverSettings object containing the wanted parameters for iterative solving
        """
        self.number_of_successful_iterations = 0
        self.number_of_unsuccessful_iterations = 0
        self.settings = settings
        super().__init__(linear_system)

    def __str__(self):
        return f"{super().__str__()}\n   Iterative"
