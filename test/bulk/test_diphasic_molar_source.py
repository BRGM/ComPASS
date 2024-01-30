import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import ComPASS
from ComPASS.utils.units import *
import MeshTools as MT
import ComPASS.io.mesh as io
import sys
from ComPASS.postprocess import postprocess
from ComPASS.timeloops import standard_loop

# ******************************************************************************
# Setting up initial values

omega_matrix = 0.15  # reservoir porosity
k_matrix = 1e-12  # reservoir permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K
Lx, Ly, Lz = 100.0, 5.0, 5.0
Ox, Oy, Oz = 0, -Ly / 2.0, -Lz / 2.0
nx, ny, nz = 100, 5, 5
p0 = 1.0e6
T0 = 300.0
molar_flux = 5e-2

# ******************************************************************************
# Grid
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

# ******************************************************************************
# The simulation object

simulation = ComPASS.load_physics("diphasic")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)


# ******************************************************************************
def select_dirichlet():
    xyz = simulation.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    return (x == x.min()) | (x == x.max())


def cell_molar_sources():
    centers = simulation.compute_global_cell_centers()
    res = np.zeros((len(centers), 2))
    res[:, 1] = molar_flux
    return res


def cell_heat_source():
    # energy inflow is approximated using p0, T0
    cell_heat_flux = molar_flux * simulation.liquid_molar_enthalpy(p0, T0, [0.0, 1.0])
    centers = simulation.compute_global_cell_centers()
    res = np.zeros(len(centers))
    res[:] = cell_heat_flux
    return res


# ******************************************************************************
simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet,
    cell_porosity=omega_matrix,
    cell_permeability=k_matrix,
    cell_thermal_conductivity=K,
    # don't forget to set cell_heat_source when setting cell_molar_sources
    cell_molar_sources=cell_molar_sources,
    cell_heat_source=cell_heat_source,
)

# copy is important, it is a pointer
final_molar_sources = simulation.all_molar_sources_vol().copy()
# reset the molar sources to begin without sources
molar_source = simulation.all_molar_sources_vol()
molar_source *= 0.0

# ******************************************************************************
initial_state = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)

for states in (
    simulation.dirichlet_node_states(),
    simulation.all_states(),
):
    states.set(initial_state)

# *************************************************************************************
init_dt = 0.01 * hour
final_time = 1 * hour


def increase_molar_sources(tick):
    tick_mol_sources = tick.time / final_time * final_molar_sources
    molar_source = simulation.all_molar_sources_vol()
    # to modify the simulation.all_molar_sources_vol() object, use += or *=
    molar_source *= 0
    molar_source += tick_mol_sources


standard_loop(
    simulation,
    initial_timestep=init_dt,
    final_time=final_time,
    output_period=10 * minute,
    iteration_callbacks=[increase_molar_sources],
)

# # For visualization, results can be postprocessed here
# simulation.postprocess(
#     convert_temperature=True,
# )

abs_centers = abs(simulation.compute_cell_centers())
where = np.logical_and(abs_centers[:, 1] < 0.1, abs_centers[:, 2] < 0.1)
x_cell_pressure = simulation.cell_states().p[where]

plt.plot(x_cell_pressure)
plt.xlabel("X cells")
plt.ylabel("Pressure (Pa)")
plt.title("Overpressure profile along X cells")
plt.savefig(
    "cells_overpressure.png",
    format="png",
)
