#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import MeshTools as MT
import MeshTools.GridTools as GT
from MeshTools.GridTools import GridInfo as Grid

# We must load mpi wrapper first so that MPI is  initialized before calling PETSC_Initialize
from . import mpi

from . import runtime
from .utils import filenames as utils_filenames
from . import dumps
from . import messages
from . import newton
from . import timestep

from .RawMesh import RawMesh
from .runtime import to_output_directory, set_output_directory_and_logfile

from ._kernel import get_kernel, load_eos, common_wrapper


# __variables__ with double-quoted values will be available in setup.py:
__version__ = "4.1.0"


def __getattr__(name):
    try:
        return getattr(common_wrapper, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
