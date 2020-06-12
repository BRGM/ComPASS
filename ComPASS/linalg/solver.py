#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple

IterativeSolverSettings = namedtuple(
    "IterativeSolverSettings",
    ["tolerance", "max_iterations", "gmres_restart"],
    defaults=[None],
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
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.linear_system = linear_system


class DirectSolver(LinearSolver):
    def __init__(self, linear_system):
        super().__init__(linear_system)
        self.activate_direct_solver = True


class IterativeSolver(LinearSolver):
    def __init__(self, linear_system, settings):
        """
        :param settings: An IterativeSolverSettings object containing the wanted parameters for iterative solving
        """
        super().__init__(linear_system)
        self.settings = settings
        self.activate_direct_solver = False
