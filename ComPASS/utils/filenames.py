#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
from .. import mpi

default_case_name = "compass"


def create_directories(path):
    if mpi.is_on_master_proc:
        if not os.path.exists(path):
            os.makedirs(path)
        assert os.path.isdir(path)


def output_directory(case_name=None, rootname=None):
    if rootname is None:
        rootname = os.getcwd()
    if case_name is None:
        case_name = default_case_name
    else:
        case_name = os.path.splitext(os.path.basename(case_name))[0]
    output = os.path.join("output-" + os.path.splitext(os.path.basename(case_name))[0])
    output = os.path.abspath(output)
    # master proc manages directory creation
    create_directories(output)
    return output


def output_directory_and_logfile(case_name=None):
    output = output_directory(case_name)
    if case_name is None:
        case_name = default_case_name
    case_name = os.path.splitext(os.path.basename(case_name))[0]
    logfile = os.path.join(output, case_name + ".log")
    return output, logfile
