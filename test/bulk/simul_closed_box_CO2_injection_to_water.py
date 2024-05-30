# This script is a test case for the flash routine put in place
# to simulate two phases (liquid, gas) two component (CO2, H2O) flows in porous media.
# It represents a simple geometry corresponding to a 1-D flow filled with liquid H2O.
# There is a CO2 injection (debit_vol) represented through Neumann
#
#                   |--------------------------|
#  debit_vol  CO2  =>       P_init water       |
#                   |--------------------------|
#
# Adding quantity to the model conduct to the gas phase appearance (there is no water vapour,
# so only CO2 in the gas phase); together with the liquid phase (diphasic_context)
# gas_context = 1; liquid_context = 2; diphasic_context = 3
#%%
import ComPASS
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

H = 10
dh = 0.1
nxy = 1

T_init = 280
P_init = 10.0 * bar
debit_vol = 1e-8  # (m/s)
Composition_inj = [1, 0]  # inject CO2, no water

omega = 0.15  # reservoir porosity
kh = 5e-14  # reservoir horizontal permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K

# Load EoS
s = simulation = ComPASS.load_physics("diphasicCO2")

ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0.0)

# Load grid
nh = int(H // dh)
dh = H / nh
L = dh * nxy / 2
grid = ComPASS.Grid(
    shape=(nxy, nxy, nh),
    extent=(2 * L, 2 * L, H),
    origin=(-L, -L, 0),
)

# energie init
E_inj = s.cpp_gas_molar_enthalpy(P_init, T_init, Composition_inj)
dens_inj = s.cpp_gas_molar_density(P_init, T_init, Composition_inj)
# print(dens_inj, E_inj)


def permeability():
    return kh * np.eye(3, dtype="d")


simulation.init(
    mesh=grid,
    cell_porosity=omega,
    cell_permeability=permeability,
    cell_thermal_conductivity=K,
)

init_state = simulation.build_state(simulation.Context.liquid, P_init, T_init)


def set_states(states, X):
    states.set(X)


set_states(simulation.node_states(), init_state)
set_states(simulation.cell_states(), init_state)

# simulation.reset_dirichlet_nodes(on_zmax(grid))
# set_states(simulation.dirichlet_node_states(), init_state)

Neumann = ComPASS.NeumannBC(
    [
        dens_inj * debit_vol * Composition_inj[0],
        dens_inj * debit_vol * Composition_inj[1],
    ],
    dens_inj * E_inj * debit_vol,
)
face_centers = simulation.face_centers()
simulation.set_Neumann_faces(on_zmin(grid)(face_centers), Neumann)


# *************************************************************************************

# output a VTU mesh file with point properties
def write_properties(filename):
    petrophysics = simulation.petrophysics()
    # pointdata – a dictionnary of point based properties
    pointdata = {
        "dirichlet": simulation.dirichlet_nodes(),
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": simulation.temperature_dirichlet_values(),
    }
    # celldata – a dictionnary of cell based properties
    celldata = {
        "Pressure": simulation.cell_states().p,
        "phi": petrophysics.cell_porosity,
    }
    celldata.update(
        tensor_coordinates(petrophysics.cell_permeability, "k", diagonal_only=True)
    )
    celldata.update(
        tensor_coordinates(
            petrophysics.cell_thermal_conductivity, "K", diagonal_only=True
        )
    )
    io.write_mesh(simulation, filename, pointdata=pointdata, celldata=celldata)


# write_properties("simulation_mesh2")


def is_correct(*args):
    x = simulation.node_states()
    correct = np.all(x.S[x.context == s.Context.liquid.value, 0] == 0)
    assert correct


# lsolver = linear_solver(simulation, direct=False)
newton = Newton(simulation, 1e-6, 20, linear_solver(simulation))

Nb_snap = 22

final_time = 100 * day
simulation.standard_loop(
    time_step_manager=TimeStepManager(
        initial_timestep=1 * day,
        increase_factor=1.3,
        decrease_factor=0.5,
        minimum_timestep=10,
    ),
    newton=newton,
    final_time=final_time,
    nb_output=Nb_snap,
    newton_iteration_callbacks=[is_correct],
)

# Simulation results can be directly postprocessed here
# simulation.postprocess()
