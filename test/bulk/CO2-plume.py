# This script is a test case for the flash routine put in place
# to simulate two phases (liquid, gas) One component (CO2) flows in porous media.
# It represents a simple geometry corresponding to a 2D plume
# The gas phase is injecte on the bottom left part of the domain with co2_flux (m/s)
# An optional caprock has thickness hcap
# The reservoir has thickness hres
#                      <----------- L ------------>
#  -                   |--------------------------|
#  |                   |                          |
# hcap                 |       caprock (kcap)     |
#  |                   |                          |
#  -                   |--------------------------|
#  |                   |                          |
#  |                   |                          |
# hres                 |      reservoir (kres)    |
#  |                   |                          |
#  |                   |                          |
#  |   -    co2_flux=> |                          |
#  |  hinj  co2_flux=> |                          |
#  -   -    co2_flux=> |--------------------------|
#
# %%

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.grid import on_zmin
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.units import *
import numpy as np
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

from ComPASS.postprocess import postprocess
from ComPASS.timeloops import standard_loop
from ComPASS.utils.various import tensor_coordinates
import ComPASS.io.mesh as io
from ComPASS.utils.mesh import extrude_triangle_strips, extrude_quad_strips

cell_type = "hex"
assert cell_type in ["hex", "prism"]

L = 300
hres = 10  # thickness of reservoir in meters
hcap = hres  # thickness of caprock in meters
hinj = 5  # thickness of injection zone in meters (must be in reservoir)
nx = 31  # number of points (number of cells + 1) along Ox
nz = 11  # number of points (number of cells + 1) along Oz

T_init = 330  # temperature of the reservoir in K
P_init = 100.0 * bar  # pressure at the bottom of the reservoir
co2_flux = 1e-7  # (m/s)

assert hres > 0
assert hcap >= 0
H = hres + hcap
dz = H / (nz - 1)
assert 0 < dz <= hinj  # at least one injection face
assert hinj < H
assert hcap >= 0
assert hinj < (H - hcap)

omega = 0.15  # porosity
kres = 5e-14  # reservoir permeability in m^2
kcap = 1e-16  # caprock permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K

# Load EoS
simulation = ComPASS.load_physics("diphasicCO2")

ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(9.81)

mesh = None
if mpi.is_on_master_proc:
    mesh = {"prism": extrude_triangle_strips, "hex": extrude_quad_strips,}[
        cell_type
    ](L, nx, dz, nbstrips=nz - 1)


def permeability():
    z = simulation.compute_global_cell_centers()[:, 2]
    k = np.zeros(z.shape[0], dtype=np.double)
    k[:] = kcap
    k[z < (H - hcap)] = kres
    return k


simulation.init(
    mesh=mesh,
    cell_porosity=omega,
    cell_permeability=permeability,
    cell_thermal_conductivity=K,
)

if mesh is not None:
    del mesh
    mesh = None

X0 = simulation.build_state(simulation.Context.liquid, P_init, T_init)

simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)
z = simulation.vertices()[:, 2]
simulation.reset_dirichlet_nodes(z < 0.5 * dz)

newton = Newton(simulation, 1e-6, 20, linear_solver(simulation))

mpi.master_print("\n" + "=" * 100 + "\n")
mpi.master_print(" " * 40 + "hydrostatic equilibrium")
mpi.master_print("\n" + "=" * 100 + "\n")

final_time = 100 * year
simulation.standard_loop(
    time_step_manager=TimeStepManager(
        initial_timestep=1 * day,
        increase_factor=4,
        decrease_factor=0.5,
        minimum_timestep=1,
    ),
    newton=newton,
    final_time=final_time,
    no_output=True,
)

mpi.master_print("\n" + "=" * 100 + "\n")
mpi.master_print(" " * 40 + "injection")
mpi.master_print("\n" + "=" * 100 + "\n")

x = simulation.vertices()[:, 0]
simulation.reset_dirichlet_nodes(x >= L)

# rescale injection flux such that it is consistent with discretization
nb_inj_faces = int(hinj / dz + 0.5)
co2_flux *= hinj / (nb_inj_faces * dz)
Composition_inj = [1, 0]  # inject CO2, no water
# energie init CO2
E_inj = simulation.cpp_gas_molar_enthalpy(P_init, T_init, Composition_inj)
dens_inj = simulation.cpp_gas_molar_density(P_init, T_init, Composition_inj)
face_centers = simulation.face_centers()
fx = face_centers[:, 0]
fz = face_centers[:, 2]
inj_faces = (fx <= 0) & (fz <= hinj)
Neumann = ComPASS.NeumannBC(
    [
        dens_inj * co2_flux * Composition_inj[0],
        dens_inj * co2_flux * Composition_inj[1],
    ],
    dens_inj * E_inj * co2_flux,
)
simulation.set_Neumann_faces(inj_faces, Neumann)

final_time = 2 * year
simulation.standard_loop(
    time_step_manager=TimeStepManager(
        initial_timestep=1 * hour,
        increase_factor=1.2,
        decrease_factor=0.8,
        minimum_timestep=1,
    ),
    newton=newton,
    final_time=final_time,
    output_period=final_time / 10,
)

# Simulation results can be directly postprocessed here
simulation.postprocess()
