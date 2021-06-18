import numpy as np

import ComPASS
from ComPASS.utils.units import *

from my_kr import kr_functions

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")

ztop = 2000
zinterface = 1000
zbottom = 0
# FIXME: for a weird reason timestep collapses around 37.428y with nb_layers=20
nb_layers = 40  # number of horizontal layers
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

face_centers = simulation.face_centers()
bottom_face = face_centers[:, 2] <= 0
bottom_nodes = simulation.facenodes(face_centers[:, 2] <= 0)
node_states = simulation.node_states()
hg = simulation.gas_molar_enthalpy
hl = simulation.liquid_molar_enthalpy
rhog = simulation.gas_molar_density
rhol = simulation.liquid_molar_density


def update_flux(*args):
    if len(bottom_nodes) == 0:
        return  # Nothing to do, may happen in parallel
    p = node_states.p[bottom_nodes]
    T = node_states.T[bottom_nodes]
    Sg = node_states.S[bottom_nodes, 0]
    rhogSg = rhog(p, T) * Sg
    gmf = rhogSg / (rhogSg + rhol(p, T) * (1 - Sg))  # gmf = gas mass fraction
    Neumann.heat_flux = -mass_flux * np.mean(gmf * hg(p, T) + (1 - gmf) * hl(p, T))
    simulation.clear_all_neumann_contributions()
    simulation.set_Neumann_faces(bottom_face, Neumann)


update_flux()

simulation.standard_loop(
    initial_timestep=day,
    output_period=30 * day,
    final_time=40 * year,
    iteration_callbacks=[update_flux,],
)

simulation.postprocess()
