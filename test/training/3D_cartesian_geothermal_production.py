# Documentation : https://compass.gitlab.io/v4/doc/
import numpy as np
import ComPASS
from ComPASS.utils.units import *  # contains MPa, degC2K, km, year...
from ComPASS.utils.grid import grid_center

# reservoir geometry
Lx, Ly, Lz = 2000.0, 2000.0, 500.0
Ox, Oy, Oz = -1000.0, -1000.0, -1500.0
nx, ny, nz = 21, 21, 20  # nz must be even


# reservoir petrophysics
k_reservoir = 1e-13  # reservoir permeability in m^2
omega_reservoir = 0.15  # reservoir porosity
K_reservoir = 2  # reservoir bulk thermal conductivity in W/m/K
# well
Qm = 300.0 * ton / hour  # well production flowrate
# fracture
k_fracture = 1e-12  # fracture permeability in m^2
omega_fracture = 0.5  # fracture porosity
K_fracture = 2  # fracture bulk thermal conductivity in W/m/K


# -------------------------------------------------------------------
# load the water2ph physics : it contains the water component
# which can be in liquid and/or gas phase
simulation = ComPASS.load_physics("water2ph")


# -------------------------------------------------------------------
# create a Cartesian grid
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)


# -------------------------------------------------------------------
# define a vertical well at grid center
def wells_factory(grid):
    def make_well():
        Cx, Cy, _ = grid_center(grid)
        well = simulation.create_vertical_well((Cx, Cy))  # define the geometry
        well.operate_on_flowrate = (
            Qm,
            1.0 * bar,
        )  # target flow rate in kg/s, threshold pressure
        well.produce()  # define the role of the well
        return (well,)  # must be a list

    return make_well


# -------------------------------------------------------------------
# Fracture factory
# set the fracture thickness (global variable)
fracture_thickness = 1  # m
simulation.set_fracture_thickness(fracture_thickness)
# select horizontal fault in the middle of the simulation domain
def fractures_factory(grid):
    def select_fractures():
        # this selection needs nz to be an even number
        face_centers = simulation.compute_global_face_centers()
        _, _, Cz = grid_center(grid)
        dz = grid.extent[2] / grid.shape[2]
        return np.abs(face_centers[:, 2] - Cz) < 0.25 * dz

    return select_fractures


# -------------------------------------------------------------------
# Initialize the regionalized values and distribute the domain
simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    wells=wells_factory(grid),
    fracture_faces=fractures_factory(grid),
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    fracture_porosity=omega_fracture,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=K_fracture,
)


# -------------------------------------------------------------------
# Initialize the domain with liquid phase
# at p_reservoir pressure and T_reservoir temperature.
top_pressure = 20.0 * MPa
T_reservoir = degC2K(70.0)
# First construct the state (at thermodynamic equilibrium)
X0 = simulation.build_state(simulation.Context.liquid, p=top_pressure, T=T_reservoir)
# then apply this state everywhere to init the domain (constant values)
simulation.all_states().set(X0)
# also apply this state as Dirichlet values
simulation.dirichlet_node_states().set(X0)


# modify the initial and Dirichlet values to apply hydrostatic pressure
bottom_pressure = top_pressure + simulation.get_gravity() * 900.0 * grid.extent[2]


def set_linear_gradients_state(states, z):
    def linear_gradients(bottom_value, top_value, depth, z):
        return top_value - ((bottom_value - top_value) / depth) * z

    states.p[:] = linear_gradients(bottom_pressure, top_pressure, grid.extent[2], z)


set_linear_gradients_state(simulation.all_states(), simulation.all_positions()[:, 2])
set_linear_gradients_state(
    simulation.dirichlet_node_states(), simulation.vertices()[:, 2]
)


# -------------------------------------------------------------------
# execute the time loop
simulation.standard_loop(
    initial_timestep=day,
    final_time=20 * year,
    output_period=year,
)


# -------------------------------------------------------------------
# some postprocesses, it allows to visualize with Paraview
simulation.postprocess()
