#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

try:
    from scipy.integrate import ode
except ImportError:
    print("Install scipy for ode or implement a simple numerical integration routine.")
    raise

import ComPASS
from ComPASS.utils.units import *

ComPASS.load_physics("water2ph")

H = 3 * km
depth_boiling = 1.5 * km

depth = np.linspace(0, depth_boiling, 200)
g = ComPASS.get_gravity()
dpdz = lambda _, p: ComPASS.liquid_volumetric_mass_density(p, ComPASS.Tsat(p)) * g
p = [1 * bar]
r = ode(dpdz).set_integrator("lsoda")
r.set_initial_value(p[0], depth[0])
p.extend([float(r.integrate(di)) for di in depth[1:]])
p = np.array(p)
Tsat = ComPASS.Tsat(p)
pdb, Tdb = p[-1], Tsat[-1]
print(
    "boiling point-depth conditions at ",
    depth[-1],
    "m:",
    pdb / MPa,
    "MPa",
    K2degC(Tdb),
    "degC",
)
rho = simulation.liquid_volumetric_mass_density(pdb, Tdb)
print(
    "bottom hydrostatic conditions",
    (pdb + rho * g * (H - depth[-1])) / MPa,
    "MPa",
    K2degC(Tdb),
    "degC",
)
