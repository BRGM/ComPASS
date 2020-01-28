#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.wells import create_vertical_well
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.io.mesh as io


H    = 1000
df   = 0.5
rw   = 0.1
epsilon = 0.0001
pres = 5. * MPa                  # initial reservoir pressure
Tres = degC2K( 70. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 30. )        # injection temperature - convert Celsius to Kelvin degrees
omega_reservoir = 0.15            # reservoir porosity
k_fracture  = 1E-11 
k_reservoir = 1E-14               # reservoir permeability in m^2
K_reservoir = 2                   # bulk thermal conductivity in W/m/K
qw = 0.1
mu = 1.E-3
rho = 1.E3

Lx = Ly = Lz = 2*H
Ox = Oy = Oz = -H
n = 20

ComPASS.set_output_directory_and_logfile(f"cemracs_4_1_{n}")

simulation = ComPASS.load_eos('water2ph')
simulation.set_gravity(0)
simulation.set_fracture_thickness(df)

grid = ComPASS.Grid(
    shape = (n, n, n/10),
    extent = (Lx, Ly, Lz),
    origin = (Ox, Oy, Oz),
)

def make_well():
    well = create_vertical_well(simulation, (0, 0), rw)
    well.operate_on_flowrate = (2*H+k_fracture/k_reservoir*df)*qw , 0. * MPa
    # well.operate_on_pressure = 1. * MPa, qmax
    well.produce()
    return [well]

def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    x, y = vertices[:,0], vertices[:,1]
    return (x <= -H + epsilon) | (x >= H - epsilon) | (y <= -H + epsilon) | (y >= H - epsilon)

def select_fracture():
    face_centers = simulation.compute_global_face_centers()
    z = face_centers[:,2]
    return np.abs(z)<epsilon

simulation.init(
    mesh = grid,
    set_dirichlet_nodes = select_dirichlet_nodes,
    wells = make_well,
    fracture_faces = select_fracture,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = K_reservoir,
    fracture_porosity = omega_reservoir,
    fracture_permeability = k_fracture,
    fracture_thermal_conductivity = 2.
)

from solution_test1 import Solution

sol = Solution(pres, qw, mu, k_reservoir, rho, rw)

# -- Set initial state and boundary conditions
initial_state = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state) # will init all variables: context, states...
# vertices is an array with shape (nb_vertices, 3)
vertices = simulation.vertices()
x, y = vertices[:,0], vertices[:,1]
dirichlet.p[:] = sol(x, y) # override the pressure value with the analytical solution

# output the initial state before the simulation is run
io.write_mesh(
    simulation,
    f"cemracs_4_1_{n}-initial_state",
    pointdata={
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": simulation.temperature_dirichlet_values(),
        "initial pressure": simulation.node_states().p,
        "initial temperature": simulation.node_states().T,
    },
    celldata={
        "initial pressure": simulation.cell_states().p,
        "initial temperature": simulation.cell_states().T,
    },
)

standard_loop(
    simulation,
    initial_timestep = 1 , final_time = year ,
    output_period = year,
)

#Compute L2 discrete error at the center
centers  = simulation.cell_centers()
usol_c   = simulation.cell_states().p
centers_x, centers_y = centers[:,0], centers[:,1]
ua_sol=sol(centers_x,centers_y)
error_c= np.linalg.norm(usol_c-ua_sol)/np.linalg.norm(ua_sol)
print("-------------------------------------------------------------------")
print("cell_centers:\t"+ "L2D-error: ", error_c, ", DOFs: ", centers.shape[0])
