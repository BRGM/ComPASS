#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
from ComPASS.timestep_management import FixedTimeStep

p0 = 0. * bar              # dummy pressure no gravity
Tmean = 10.                # average surface temperature
deltaT = 10.               # seasonal amplitude
bottom_heat_flux = 0.08    # W/m2                                  
k_matrix = 1E-18           # permeability - not relevant here but cannot be 0
phi_matrix = 0.15          # porosity - not really relevant either as long as
                           # fluid and rock volumetric heat capacities are the same
K_matrix = 2.              # bulk thermal conductivity in W/m/K
rhor = 2200.               # rock density kg/m^3
cpr = 800.                 # rock pecific heat capacity J/K/kg

H = 100.                   # column height
nx, ny, nz = 1, 1, 200     # discretization
dz = H/nz

ComPASS.activate_direct_solver = True
ComPASS.load_eos('linear_water')
ComPASS.set_rock_volumetric_heat_capacity(rhor*cpr)
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.volumetric_heat_capacity = rhor*cpr
ComPASS.set_gravity(0)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (dz, dz, H),
    origin = (-0.5*dz, -0.5*dz, -H),
)

def top_nodes():
    return ComPASS.global_vertices()[:, 2] >= 0
    
def set_node_flags():
    ComPASS.global_nodeflags()[:] = np.asarray(top_nodes(), dtype=np.int)

ComPASS.init(
    mesh = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    cell_thermal_conductivity = K_matrix,
    set_dirichlet_nodes = top_nodes,
    set_global_flags = set_node_flags,
)

def set_initial_states(states, z):
    states.context[:] = 1
    states.p[:] = p0
    states.T[:] = Tmean - ( bottom_heat_flux / K_matrix ) * z
    states.S[:] = 1.
    states.C[:] = 1.
set_initial_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,2])
set_initial_states(ComPASS.node_states(), ComPASS.vertices()[:,2])
set_initial_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,2])

def set_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux
    face_centers = ComPASS.face_centers()   
    ComPASS.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann) 
set_boundary_heat_flux()

# locate dirichlet nodes - not mandatory
# we could have identified different regions using nodeflags
dirichlet_nodes = np.nonzero(ComPASS.nodeflags())[0]
dirichlet_T = ComPASS.dirichlet_node_states().T

def change_surface_temperature(n, t):
    dirichlet_T[dirichlet_nodes] = Tmean + deltaT * np.sin(t * (2*np.pi / year))

final_time = 5 * year
output_period = year / 12
standard_loop(
   final_time = final_time,
   time_step_manager = FixedTimeStep(5. * day),
   output_period = output_period,
# iteration callbacks are function of the form f(n, t) 
# where n is the iteration number and t is time
# they are called at the end of each successful iteration
# you can put as many of them
   iteration_callbacks = [change_surface_temperature],
)
