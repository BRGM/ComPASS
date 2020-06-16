#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


class LinearSolverFailure(Exception):
    def __init__(self, reason):
        self.reason = reason


class DirectSolverFailure(LinearSolverFailure):
    def __str__(self):
        return f"Direct solver failed with reason :\n{self.reason}"


class IterativeSolverFailure(LinearSolverFailure):
    def __init__(self, reason, nit):
        super().__init__(reason)
        self.nit = nit

    def __str__(self):
        return f"Iterative solver failed after {self.nit} iterations with reason: {self.reason}"
