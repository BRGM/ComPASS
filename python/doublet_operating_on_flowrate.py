#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

grid = ComPASS.Grid(
    shape = (31, 21, 1),
    extent = (3000., 2000., 20),
    origin = (-1500, -1000, -1600),
)

def make_wells():
	producer = doublet_utils.make_well((-500, 0))
	producer.operate_on_flowrate = 300 * ton / hour, 1E5
	producer.produce()
	injector = doublet_utils.make_well((500, 0))
	injector.operate_on_flowrate = 300 * ton / hour, 30E6
	injector.inject(degC2K(30))
	return [producer, injector]

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    wells = make_wells
)

standard_loop(final_time = 30 * year, output_period = year)
