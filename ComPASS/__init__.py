#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# We must load mpi wrapper first so that MPI is  initialized before calling PETSC_Initialize
from . import mpi

from MeshTools.GridTools import GridInfo as Grid
from .RawMesh import RawMesh
from .runtime import to_output_directory, set_output_directory_and_logfile
from ._kernel import get_kernel, load_eos, common_wrapper

import verstr

try:
    from . import _version

    __version__ = verstr.verstr(_version.version)
except ImportError:
    __version__ = None


def __getattr__(name):
    try:
        return getattr(common_wrapper, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
