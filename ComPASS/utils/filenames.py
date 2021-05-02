#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
from pathlib import Path
from .. import mpi

default_case_name = "compass"


def create_directories(path):
    if mpi.is_on_master_proc:
        if not os.path.exists(path):
            os.makedirs(path)
        assert os.path.isdir(path)


def output_directory(case_name=None, rootname=None, process_case_name=True):
    if rootname is None:
        rootname = os.getcwd()
    if case_name is None:
        case_name = default_case_name
    else:
        if process_case_name:
            case_name = os.path.splitext(os.path.basename(case_name))[0]
        output = os.path.join("output-" + case_name)
        output = os.path.abspath(output)
    # master proc manages directory creation
    create_directories(output)
    return output


def output_version_info(path):
    if mpi.is_on_master_proc:  # output directory is created by master proc
        import platform

        root = Path(__file__).parent.parent
        with open(os.path.join(path, "version_info"), "w") as f:
            with (root / "version_info").open() as build_info:
                for line in build_info:
                    f.write(line)
            print(
                f"running python {platform.python_version()} on {platform.node()}",
                file=f,
            )
            print(f"package location: {str(root)}", file=f)
            print(f"system info: {platform.platform()}", file=f)


def output_directory_and_logfile(
    case_name=None, process_case_name=True, add_version_info=True
):
    output = output_directory(case_name)
    if add_version_info:
        output_version_info(output)
    if case_name is None:
        case_name = default_case_name
    elif process_case_name:
        case_name = os.path.splitext(os.path.basename(case_name))[0]
    logfile = os.path.join(output, case_name + ".log")
    return output, logfile
