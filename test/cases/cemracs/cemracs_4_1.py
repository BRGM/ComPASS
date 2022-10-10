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
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timeloops import standard_loop
import ComPASS.io.mesh as io
import ComPASS.mpi as mpi
from ComPASS.mpi import MPI  # underlying mpi4py.MPI
from ComPASS.properties.utils import constant_physical_property
from ComPASS.properties.densities import build_pure_phase_volumetric_mass_density


H = 1000
df = 0.5
rw = 0.1
epsilon = 0.0001
production = True
if production:
    pres = 5.0 * MPa  # initial reservoir pressure, MPa=10^6
else:  # injection
    pres = 20.0 * MPa  # initial reservoir pressure, MPa=10^6
    Tinjection = degC2K(
        30.0
    )  # injection temperature - convert Celsius to Kelvin degrees
Tres = degC2K(
    70.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
omega_reservoir = 0.15  # reservoir porosity
k_fracture = 1e-11
kres = 1e-14  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K
qw = 1.0
mu = 1.0e-3
rho = 1.0e3

Lx = Ly = Lz = 2 * H
Ox = Oy = Oz = -H
n = 10

ComPASS.set_output_directory_and_logfile(f"cemracs_4_1_{n}")

simulation = ComPASS.load_physics("linear_water")
simulation.set_molar_density_functions(
    build_pure_phase_volumetric_mass_density(
        specific_mass=rho,
        compressibility=0,
        thermal_expansivity=0,
    ),
)
simulation.set_viscosity_functions(constant_physical_property(mu))

# fluid_properties.volumetric_heat_capacity = rhor*cpr # not relevant here - use default value
# simulation.set_rock_volumetric_heat_capacity(rhor*cpr) # not relevant here - use default value

simulation.set_gravity(0)
simulation.set_fracture_thickness(df)

grid = ComPASS.Grid(
    shape=(n, n, n),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)


def make_well():
    well = simulation.create_vertical_well((0, 0), rw)
    if production:
        well.operate_on_flowrate = (Lz + (k_fracture / kres) * df) * qw, 0.0 * MPa
        well.produce()
    else:  # injection
        well.operate_on_flowrate = (
            (Lz + (k_fracture / kres) * df) * qw,
            pres + 100.0 * MPa,
        )
        well.inject(Tinjection)
    return [well]


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    x, y = vertices[:, 0], vertices[:, 1]
    return (
        (x <= -H + epsilon)
        | (x >= H - epsilon)
        | (y <= -H + epsilon)
        | (y >= H - epsilon)
    )


def select_fracture():
    face_centers = simulation.compute_global_face_centers()
    z = face_centers[:, 2]
    return np.abs(z) < epsilon


def build_anisotropic_permeability():
    centers = simulation.compute_global_cell_centers()
    result = np.zeros((centers.shape[0], 3, 3), dtype=np.double)
    result[..., 0, 0] = kres
    result[..., 1, 1] = kres
    return result


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    wells=make_well,
    fracture_faces=select_fracture,
    cell_porosity=omega_reservoir,
    cell_permeability=build_anisotropic_permeability,
    cell_thermal_conductivity=K_reservoir,
    fracture_porosity=omega_reservoir,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=2.0,
)

from solution_test1 import Solution

sol = Solution(pres, qw, mu, kres, rho, rw)

# -- Set initial state and boundary conditions
initial_state = simulation.build_state(p=pres, T=Tres)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)  # will init all variables: context, states...
# vertices is an array with shape (nb_vertices, 3)
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
dirichlet.p[:] = sol(x, y)  # override the pressure value with the analytical solution

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

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

standard_loop(
    simulation,
    newton=newton,
    initial_timestep=0.1,
    final_time=40,
    output_period=year,  # year ,
)

rank = mpi.proc_rank
ncells = simulation.nb_cells_own()[rank]
centers = simulation.cell_centers()[:ncells]
usol_c = simulation.cell_states().p[:ncells]
centers_x, centers_y = centers[:, 0], centers[:, 1]
ua_sol = sol(centers_x, centers_y)
sum_squares = lambda a: np.sum(a**2)
error_c = np.array([sum_squares(usol_c - ua_sol), sum_squares(ua_sol)], np.double)
if not mpi.is_on_master_proc:
    MPI.COMM_WORLD.Reduce(
        [error_c, MPI.DOUBLE], None, op=MPI.SUM, root=mpi.master_proc_rank
    )
else:
    result = np.zeros_like(error_c)
    MPI.COMM_WORLD.Reduce(
        [error_c, MPI.DOUBLE],
        [result, MPI.DOUBLE],
        op=MPI.SUM,
        root=mpi.master_proc_rank,
    )
    print("-------------------------------------------------------------------")
    relative_error = np.sqrt(result[0]) / np.sqrt(result[1])
    print(f"cell_centers:\t L2D-error: {relative_error} DOFs: {np.prod(grid.shape)}")
