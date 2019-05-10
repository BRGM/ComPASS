
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager


def test__extrems__no_gravity(
    T_injection_degC = 33.0,
    flow_velocity_m_s = 1.E-6
):

    T_injection = degC2K( T_injection_degC )

    p0 = 1. * bar
    T0_degC = 5.
    T0 = degC2K( T0_degC )
    p_outlet = p0
    permeability = 1E-12 # m^2
    porosity = 0.2
    final_time = 1E8 # s

    rhow = 1000      # kg/m3
    Cw = 4200.       # J/K/kg
    Kw = 0.6         # W/m/K
    rhor = 2600.     # kg/m3
    Cr = 800.        # J/K/kg
    Kr = 2.          # W/m/K

    Keq = (1-porosity)*Kr + porosity*Kw

    # grid specifications
    extent = Lx, Ly, Lz = 1000, 10, 10
    shape = nx, ny, nz = 100, 1, 1
    origin = (0, -0.5*Ly, -0.5*Lz)

    ComPASS.load_eos('linear_water')
    ComPASS.set_gravity(0)
    ComPASS.set_output_directory_and_logfile(__file__)
    ComPASS.set_rock_volumetric_heat_capacity(rhor*Cr)

    p = ComPASS.get_fluid_properties()
    p.compressibility = 1e-8

    grid = ComPASS.Grid(shape = shape, extent = extent, origin = origin)

    def outlet_nodes():
       return ComPASS.global_vertices()[:, 0] >= Lx

    ComPASS.init(
       grid = grid,
       cell_permeability = permeability,
       cell_porosity = porosity,
       cell_thermal_conductivity = Keq,
       set_dirichlet_nodes = outlet_nodes,
    )

    def set_initial_states(states):
       states.context[:] = 1 # liquid
       states.p[:] = p0
       states.T[:] = T0
       states.S[:] = 1
       states.C[:] = 1.     # component fraction... here only one component
    for states in [ComPASS.dirichlet_node_states(),
                  ComPASS.node_states(),
                  ComPASS.cell_states()]:
       set_initial_states(states)

    def set_boundary_flux():
        Neumann = ComPASS.NeumannBC()
        specific_massflux = flow_velocity_m_s * ComPASS.molar_density(p0, T_injection) / (Ly * Lz)
        Neumann.molar_flux[:] = specific_massflux
        # energy inflow is approximated using p0
        Neumann.heat_flux = specific_massflux * ComPASS.molar_enthalpy(p0, T_injection)
        ComPASS.set_Neumann_faces(ComPASS.face_centers()[:, 0] <= 0, Neumann)
    set_boundary_flux()

    output_period = 0.1 * final_time

    # On teste a chaque pas de temps si les valeurs extremes sont bien sur les bords

    (boundary_idx,) = np.where(
        (ComPASS.vertices()[:, 0] >= Lx)
        | (ComPASS.vertices()[:, 0] <= 0)
    )

    def valid(arr, atol=1e-3):
        arr_boundary = arr[boundary_idx]
        r_min = arr_boundary.min(0) - atol
        r_max = arr_boundary.max(0) + atol
        g_min, g_max = arr.min(0), arr.max(0)
        return np.all(
            (r_min <= g_min) & (g_min <= r_max)
            & (r_min <= g_max) & (g_max <= r_max)
        )

    def valid_current(iteration, time):
        X = ComPASS.node_states()
        assert valid(X.T)
        assert valid(X.p)
        assert valid(X.C)
        ##assert valid(X.S)  # NON pas la saturation car la bulle se developpe avant la sortie

    standard_loop(
        final_time = final_time,
        time_step_manager = TimeStepManager(1 * hour, 0.2 * output_period),
        output_period = output_period, output_callbacks=[valid_current],
    )

if __name__ == '__main__':
    test__extrems__no_gravity()
    #test__extrems__no_gravity(T_injection_degC=200)
    #test__extrems__no_gravity(flow_velocity_m_s=4e-5)
    print('test ok')