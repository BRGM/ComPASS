import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.utils.grid import on_zmax

from my_kr import kr_functions

ztop = 2000
zinterface = 1000
zbottom = 0
# FIXME: for a weird reason timestep collapses around 37.428y with nb_layers=20
nb_layers = 20  # number of horizontal layers
H = ztop - zbottom
dz = H / nb_layers
dx = dy = dz
omega_top = 0.15  # overburden porosity
k_top = 5e-15  # overburden permeability
omega_bottom = 0.25  # overburden porosity
k_bottom = 1e-10  # reservoir permeability in m^2
K = 1  # bulk thermal conductivity in W/m/K
rock_heat_capacity = 1e3  # J/kg/K
rock_density = 2500  # kg/m3
mass_flux = 100.0e-6  # kg/s - here 100 kg/s/km^2
Ttop = degC2K(10)
Tinterface = degC2K(290)
Tbottom = degC2K(310)
ptop = 1.013e5  # Pa
gravity = 9.81
final_time = 70 * year
output_period = year

# -----------------------------------------------------------------------------

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("water2ph")
simulation.set_rock_volumetric_heat_capacity(rock_heat_capacity * rock_density)
simulation.set_gravity(gravity)

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

simulation.set_kr_functions(kr_functions)

# initial values
Xliq = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(Xliq)


@np.vectorize
def geotherm(z):
    if z > zinterface:
        return (Ttop - Tinterface) * (z - zinterface) / (ztop - zinterface) + Tinterface
    else:
        return (Tinterface - Tbottom) * (z - zbottom) / (zinterface - zbottom) + Tbottom


hydrostatic_profile = simulation.hydrostatic_pressure_profile(
    zbottom, ztop, 3 * nb_layers, ptop, geotherm
)
states = simulation.all_states()
z = simulation.all_positions()[:, 2]
states.p[:] = hydrostatic_profile(z)
states.T[:] = geotherm(z)

# boundary conditions
simulation.reset_dirichlet_nodes(on_zmax(grid))
Neumann = ComPASS.NeumannBC(-mass_flux, compute_heat_flux=True, nz=-1)
bottom_face = simulation.face_centers()[:, 2] <= 0
simulation.set_Neumann_faces(bottom_face, Neumann)

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 20, lsolver)

simulation.standard_loop(
    initial_timestep=day,
    output_period=output_period,
    final_time=final_time,
    newton=newton,
)

simulation.postprocess()
