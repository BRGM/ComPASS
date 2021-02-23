import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import ComPASS
from ComPASS.utils.units import *
import sys
from ComPASS.postprocess import postprocess
from ComPASS.timeloops import standard_loop
from ComPASS.properties.densities import build_pure_phase_volumetric_mass_density
from ComPASS.properties.utils import constant_physical_property

#%%Setting up initial values-----------------------------------------------------------------

omega_matrix = 0.15  # reservoir porosity
k_matrix = 1e-12  # reservoir permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K
Lx, Ly, Lz = 100.0, 5.0, 5.0
Ox, Oy, Oz = 0, -Ly / 2.0, -Lz / 2.0
nx, ny, nz = 100, 5, 5
p0 = 1.0e6
T0 = 300.0
molar_flux = 1e-3

# Grid-----------------------------------------------------------------

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

# The simulation object---------------------------------------------------
simulation = ComPASS.load_physics("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)

# Dirichlet&source-----------------
def select_dirichlet():
    xyz = simulation.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    return (x == x.min()) | (x == x.max())


def cell_molar_sources():
    centers = simulation.compute_global_cell_centers()
    res = np.zeros((len(centers), 1))
    res[:] = molar_flux
    return res


def cell_heat_source():
    # energy inflow is approximated using p0, T0
    cell_heat_flux = molar_flux * simulation.molar_enthalpy(p0, T0)
    centers = simulation.compute_global_cell_centers()
    res = np.zeros(len(centers))
    res[:] = cell_heat_flux
    return res


rhof = 1e6 / 18  # specific mass in kg/m^3
muf = 1e-3  # viscosity Pa.s
simulation.set_molar_density_functions(
    build_pure_phase_volumetric_mass_density(specific_mass=rhof),
)
simulation.set_viscosity_functions(constant_physical_property(muf))

#%% Simulation -----------------------------------------------------------------

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

#%% homogeneous reservoir initial state # top Pressure # top Temperature---------------


def set_initial_states(states):
    states.context[:] = 1
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = 1.0
    states.C[:] = 1.0


for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
]:
    set_initial_states(states)

# *************************************************************************************
init_dt = 0.01 * hour
final_time = hour

standard_loop(
    simulation,
    initial_timestep=init_dt,
    final_time=final_time,
    output_period=minute,
    nitermax=300,
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
