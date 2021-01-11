#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os

from .utils.filenames import output_directory_and_logfile

# output_directory, logfile = output_directory_and_logfile()
output_directory = None
logfile = None


def to_output_directory(filename):
    assert os.path.isdir(output_directory), f"{output_directory} is not a directory"
    return os.path.abspath(os.path.join(output_directory, filename))


def set_output_directory_and_logfile(case_name, process_case_name=True):
    global output_directory, logfile
    output_directory, logfile = output_directory_and_logfile(case_name)
