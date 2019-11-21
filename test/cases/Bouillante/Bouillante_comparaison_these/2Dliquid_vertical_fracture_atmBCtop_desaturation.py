#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 1000m in depth
# with 1 vertical fracture at Lx/2.
# Homogeneous Neumann BC at both sides and at the bottom
# atm BC at the top.
# Imposed flux at the bottom of the fracture.

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.newton import Newton, LinearSolver

ComPASS.load_eos('diphasic_FreeFlowBC')
ComPASS.activate_direct_solver = True

pure_phase_molar_fraction = [[1., 0.], [0., 1.0]]
p0 = 1. * bar              # initial reservoir pressure
T0 = degC2K( 20. )         # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
qmass = 1E-1               # 
Tbottom = degC2K( 350. )   # temperature influx
k_matrix = 1E-15           # domain permeability in m^2
phi_matrix = 0.15          # domain porosity
k_fracture = 1E-12         # fracture permeability in m^2
phi_fracture = 0.3         # fracture porosity
thermal_cond = 2.          # bulk thermal conductivity in W/m/K

H = 1000.                  # domain height
nH = 50                   # discretization
nx, ny, nz = 2*nH, 1, nH
Lx, Ly, Lz = 2*H, 0.1*H, H
Topz = -H+H-0.5

ComPASS.set_output_directory_and_logfile(__file__)

# thermodynamic functions are only available once the eos is loaded
pbottom = ComPASS.get_gravity() * H * 1000.
hbottom = ComPASS.liquid_molar_enthalpy(pbottom, Tbottom, pure_phase_molar_fraction)

freeflow_flag = 30  # do not modify this number

liquid_context = ComPASS.Context.liquid
diphasic_with_liq_outflow = ComPASS.Context.diphasic_FF_liq_outflow

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
        origin = (-0.5*Lx, -0.5*Ly, -H)
    )

    def set_global_flags():
        # nodes
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[vertices[:, -1] >= Topz] = freeflow_flag
        # freeflow faces, necessary to flag them
        facecenters = ComPASS.compute_global_face_centers()
        faceflags = ComPASS.global_faceflags()
        faceflags[facecenters[:, -1] >= Topz] = freeflow_flag

    def select_fractures():
        centers = ComPASS.compute_global_face_centers()
        xc = centers[:, 0]
        zc = centers[:, -1]
        return (xc == 0) #& (zc > -0.5 * H)

if not ComPASS.mpi.is_on_master_proc:
    grid = set_global_flags = select_fractures = None

ComPASS.init(
    mesh = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    cell_thermal_conductivity = thermal_cond,
    fracture_faces = select_fractures,
    fracture_permeability = k_fracture,
    fracture_porosity = phi_fracture,
    fracture_thermal_conductivity = thermal_cond,
    set_global_flags=set_global_flags,
)

def set_initial_states(states):
    states.context[:] = liquid_context
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = [0, 1]
    states.C[:] = [[1.0, 0.0], [0.0, 1.0]]
    states.FreeFlow_phase_flowrate[:] = 0.
for states in [ComPASS.node_states(),
               ComPASS.cell_states(),
               ComPASS.fracture_states()]:
    set_initial_states(states)

def set_boundary_fluxes():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = [0., qmass]
    Neumann.heat_flux = qmass * hbottom
    face_centers = ComPASS.face_centers()   
    bottom_fracture_edges = ComPASS.find_fracture_edges(face_centers[:, -1] <= -H)
    ComPASS.set_Neumann_fracture_edges(bottom_fracture_edges, Neumann) 
set_boundary_fluxes()

def set_FreeFlow_state(state):
    node_flags = ComPASS.nodeflags()
    # top
    state.context[node_flags == freeflow_flag] = diphasic_with_liq_outflow
    state.p[node_flags == freeflow_flag] = p0
    state.T[node_flags == freeflow_flag] = T0
    state.S[node_flags == freeflow_flag] = [0., 1.]
    state.C[node_flags == freeflow_flag] = [[1., 0.], [0., 1.]]
    state.FreeFlow_phase_flowrate[node_flags == freeflow_flag] = 0.
set_FreeFlow_state(ComPASS.node_states())

# set linear solver properties
newton = Newton(1e-6, 35, LinearSolver(1e-6, 50))


final_time = 300 * year
output_period = 0.05 * final_time
timestep = TimeStepManager(
    initial_timestep=1.*day,
    minimum_timestep=1e-3,
    maximum_timestep=output_period,
    increase_factor=1.3,
    decrease_factor=0.2,
)

standard_loop(
    final_time = final_time,
    time_step_manager = timestep,
    output_period = output_period,
    newton=newton,
)