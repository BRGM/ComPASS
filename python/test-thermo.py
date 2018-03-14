#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import ComPASS
from ComPASS.utils.units import *

ComPASS.load_eos('water2ph')

# liquid specific density
T = np.linspace(20, 400)
xsi = ComPASS.liquid_molar_density(1 * bar, degC2K(T))
plt.clf()
plt.title('liquid specific mass')
plt.grid(True)
plt.plot(T, xsi)
plt.xlabel('temperature (Celsius degree)')
plt.ylabel('liquid specific mass')
plt.savefig('rhol.png')

T = np.linspace(degC2K(20), degC2K(300))
h = ComPASS.liquid_molar_enthalpy(ComPASS.Psat(T), degC2K(T))

try:
    from scipy.integrate import ode
except ImportError:
    print('Install scipy for ode or implement a simple numerical integration routine.')
    raise
else:
    depth = np.linspace(0, 3000, 200) # 0 - 3km
    g = ComPASS.get_gravity()
    dpdz = lambda _, p: ComPASS.liquid_molar_density(p, ComPASS.Tsat(p)) * g
    p = [1 * bar]
    r = ode(dpdz).set_integrator('lsoda')
    r.set_initial_value(p[0], depth[0])
    p.extend([float(r.integrate(di)) for di in depth[1:]])
    p = np.array(p)
    plt.clf()
    plt.title('boiling point-depth curve')
    plt.grid(True)
    Tsat = K2degC(ComPASS.Tsat(p))
    plt.plot(Tsat, depth)
    plt.xlabel('temperature (Celsius degree)')
    plt.ylabel('depth (m)')
    plt.ylim(depth.max(), 0) # invert axis
    plt.savefig('bpd.png')
    plt.clf()
    plt.title('boiling point-depth pressure curve')
    plt.grid(True)
    plt.plot(p / MPa, depth)
    plt.xlabel('pressure (MPa)')
    plt.ylabel('depth (m)')
    plt.ylim(depth.max(), 0) # invert axis
    plt.savefig('bpd-pressure.png')
    print('bottom boiling conditions:', p[-1]/MPa, 'MPa', Tsat[-1], 'degC')
