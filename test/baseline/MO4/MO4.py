import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *

# from ComPASS.newton import Newton
# from ComPASS.linalg.factory import linear_solver
from ComPASS.timeloops import TimeStepManager
import MeshTools as MT

# import ComPASS.mpi as mpi
# from ComPASS.mpi import MPI # underlying mpi4py.MPI

from my_kr import kr_functions

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")

ztop = 2000
zinterface = 1000
zbottom = 0
nb_layers = 20  # number of horizontal layers
H = ztop - zbottom
dz = H / nb_layers
dx = dy = dz
omega_top = 0.15  # overburden porosity
k_top = 5e-15  # overburden permeability
omega_bottom = 0.25  # overburden porosity
k_bottom = 1e-10  # reservoir permeability in m^2
K = 1  # bulk thermal conductivity in W/m/K
rock_heat_capacity = 1e3  # kJ/kg/K
rock_density = 2500  # kg/m3
mass_flux = 0.001 * 1e-4 * dx * dy  # 100 kg/km2
Ttop = degC2K(10)
Tinterface = degC2K(290)
Tbottom = degC2K(310)
ptop = 1.013e5  # Pa
gravity = 9.81
simulation.set_gravity(gravity)


@np.vectorize
def geotherm(z):
    if z > zinterface:
        return (Ttop - Tinterface) * (z - zinterface) / (ztop - zinterface) + Tinterface
    else:
        return (Tinterface - Tbottom) * (z - zbottom) / (zinterface - zbottom) + Tbottom


def make_hydrostatic_pressure(zbottom, ztop, nz, nbsteps=100):
    assert zbottom < ztop
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    rho = simulation.liquid_molar_density
    p = 0
    pressures = [p]
    for zbot, ztop in zip(z[1:], z[:-1]):
        zeta = ztop
        dz = (ztop - zbot) / nbsteps
        for _ in range(nbsteps):
            p += gravity * rho(p, geotherm(zeta)) * dz
            zeta -= dz
        pressures.append(p)
    pressures = np.asarray(pressures)
    return lambda zeta: np.interp(zeta, z[::-1], pressures[::-1])


hydrostatic_pressure = make_hydrostatic_pressure(zbottom, ztop, 3 * nb_layers)

simulation.set_rock_volumetric_heat_capacity(rock_heat_capacity * rock_density)

grid = ComPASS.Grid(
    shape=(1, 1, nb_layers), extent=(dx, dy, H), origin=(-0.5 * dx, -0.5 * dy, zbottom)
)


def set_permeability():
    centers = simulation.compute_global_cell_centers()
    z = centers[:, 2]
    nb_cells = centers.shape[0]
    k = np.zeros(nb_cells, dtype="d")
    k[z < zinterface] = k_bottom  # isotropic
    k[z >= zinterface] = k_top  # isotropic
    return k


def set_porosity():
    centers = simulation.compute_global_cell_centers()
    z = centers[:, 2]
    nb_cells = centers.shape[0]
    omega = np.zeros(nb_cells, dtype="d")
    omega[z < zinterface] = omega_bottom  # isotropic
    omega[z >= zinterface] = omega_top  # isotropic
    return omega


simulation.init(
    mesh=grid,
    cell_porosity=set_porosity,
    cell_permeability=set_permeability,
    cell_thermal_conductivity=K,
)

# Set the kr functions after initialization
simulation.set_kr_functions(kr_functions)

vertices = simulation.vertices()

# Set initial values
Xliq = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(Xliq)


def set_states(states, z):
    states.p[:] = hydrostatic_pressure(z)
    states.T[:] = geotherm(z)


set_states(simulation.node_states(), vertices[:, 2])
set_states(simulation.cell_states(), simulation.compute_cell_centers()[:, 2])

# Set boundary conditions
# dirichlet at the top
simulation.reset_dirichlet_nodes(vertices[:, 2] >= ztop)

# Neumann at the bottom
Neumann = ComPASS.NeumannBC()
Neumann.molar_flux[:] = -mass_flux
# we approximate the specific enthalpy of the bottom-most cell
pbottom = hydrostatic_pressure(zbottom)
hbottom = simulation.liquid_molar_enthalpy(pbottom, Tbottom)
Neumann.heat_flux = -mass_flux * hbottom
face_centers = simulation.face_centers()
simulation.set_Neumann_faces(face_centers[:, 2] <= 0, Neumann)
# breakpoint()
simulation.standard_loop(initial_timestep=day, output_period=5 * day, final_time=year)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess()
