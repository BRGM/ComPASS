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

ComPASS.load_eos('water2ph')

pres = 20. * MPa
Tres = degC2K( 70. ) # convert Celsius to Kelvin degrees
Tinjection = degC2K( 30. )
Qm = 300. * ton / hour

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

def make_wells():
    interwell_distance = 1 * km
    Ox, Oy = doublet_utils.center(grid)[:2]
    producer = doublet_utils.make_well((Ox - 0.5 * interwell_distance, Oy))
    producer.operate_on_flowrate = Qm , 1. * bar
    producer.produce()
    injector = doublet_utils.make_well((Ox + 0.5 * interwell_distance, Oy))
    injector.operate_on_flowrate = Qm, pres + 100. * MPa
    injector.inject(Tinjection)
    return (producer, injector)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    set_dirichlet_nodes = doublet_utils.select_boundary_factory(grid),
    wells = make_wells,
)
doublet_utils.init_states(pres, Tres)

standard_loop(final_time = 30 * year, output_period = year)
