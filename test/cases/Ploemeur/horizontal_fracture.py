import numpy as np

import MeshTools as MT

import ComPASS
from ComPASS.timeloops import standard_loop
from ComPASS.utils.units import *
from ComPASS.utils.wells import create_vertical_well

#%% mesh


def geometric(x0, xn, a, n):
    assert a >= 1 and n > 0
    dx = np.cumprod(np.hstack([1, np.tile(a, max(n - 1, 0))]))
    assert np.all(dx > 0)
    x = np.cumsum(dx)
    assert np.all(x[1:] >= x[:-1])
    return x0 + ((xn - x0) / (x[-1] - x[0])) * (x - x[0])


# grid is centered on the middle of the segment between injection and production points
# hence wells positions are (-0.5 * interwell_distance, 0, 0) and (0.5 * interwell_distance, 0, 0)
interwell_distance = 10  # meters
grid_length = 3 * interwell_distance  # meters
grid_depth = interwell_distance  # meters
grid_height = 1.5 * interwell_distance  # meters
assert 0 < interwell_distance
assert interwell_distance < grid_length
n = 6
nb_steps_x_outside_wells = n
nb_steps_y = n
nb_steps_z = n
steps_x_ratio_outside_wells = 1.3
steps_y_ratio = 1.3
steps_z_ratio = 1.5

half_interwell_distance = 0.5 * interwell_distance
xout = geometric(
    half_interwell_distance,
    0.5 * grid_length,
    steps_x_ratio_outside_wells,
    nb_steps_x_outside_wells,
)
dx = xout[1] - xout[0]
xbetween = np.linspace(
    -half_interwell_distance,
    half_interwell_distance,
    np.ceil(interwell_distance / dx) + 1,
)
x = np.hstack([(-xout[1:])[::-1], xbetween, xout[1:]])
ypos = geometric(0, 0.5 * grid_depth, steps_y_ratio, nb_steps_y)
y = np.hstack([(-ypos[1:])[::-1], ypos])
zpos = geometric(0, 0.5 * grid_height, steps_z_ratio, nb_steps_z)
dz = zpos[1] - zpos[0]
z = np.hstack([(-zpos[1:])[::-1], zpos])

mesh = MT.grid3D(steps=(x, y, z))

MT.to_vtu(mesh, "ploemeur_mesh.vtu")

#%% simulation

rhof = 1e3  # specific mass in kg/m^3
cpf = 4200  # specific heat in J/kg/K
rhofcpf = rhof * cpf  # volumetric heat capacity
muf = 1e-3  # dynamic viscosity Pa.s
pres = 0.0 * MPa  # initial reservoir pressure
Tres = degC2K(
    10.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K(15.0)  # injection temperature - convert Celsius to Kelvin degrees
Qm = (rhof * 1e-3) / minute  # production/injection mass flowrate
k_fracture = 1e-12  # fracture permeability in m^2
omega_fracture = 0.55  # fracture porosity
K_fracture = 2  # bulk thermal conductivity in W/m/K
k_matrix = 1e-20  # matrix permeability in m^2
omega_matrix = 0.15  # matrix porosity
K_matrix = 2  # bulk thermal conductivity in W/m/K
fracture_thickness = 0.005

ComPASS.load_eos("linear_water")
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.compressibility = 1e-10  # it helps...
fluid_properties.volumetric_heat_capacity = rhofcpf
fluid_properties.dynamic_viscosity = muf

ComPASS.set_gravity(0)
ComPASS.set_fracture_thickness(fracture_thickness)
ComPASS.set_output_directory_and_logfile(__file__)


def build_wells():
    producer = create_vertical_well(
        (half_interwell_distance, 0), well_radius=0.1, zmin=-1.5 * dz, zmax=1.5 * dz
    )
    producer.operate_on_flowrate = Qm, -1e99  # mass flow rate, pressure limit (minimum)
    producer.produce()
    injector = create_vertical_well(
        (-half_interwell_distance, 0), well_radius=0.1, zmin=-1.5 * dz, zmax=1.5 * dz
    )
    injector.operate_on_flowrate = -Qm, 1e99  # mass flow rate, pressure limit (maximum)
    injector.inject(Tinjection)
    return (producer, injector)


def select_dirichlet_nodes():
    vertices = ComPASS.global_vertices()
    return (vertices[:, 2] == z.min()) | (vertices[:, 2] == z.max())


def select_fractures():
    face_centers = ComPASS.compute_global_face_centers()
    return face_centers[:, 2] == 0


ComPASS.init(
    mesh=mesh,
    wells=build_wells,
    cell_porosity=omega_matrix,
    cell_permeability=k_matrix,
    cell_thermal_conductivity=K_matrix,
    fracture_faces=select_fractures,
    fracture_porosity=omega_fracture,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=K_fracture,
    set_dirichlet_nodes=select_dirichlet_nodes,
)


def set_states(states):
    states.context[:] = 1
    states.p[:] = pres
    states.T[:] = Tres
    states.S[:] = 1.0
    states.C[:] = 1.0


for states in [
    ComPASS.dirichlet_node_states(),
    ComPASS.node_states(),
    ComPASS.fracture_states(),
    ComPASS.cell_states(),
]:
    set_states(states)

final_time = 2e4
output_period = 0.05 * final_time
ComPASS.set_maximum_timestep(0.5 * output_period)
ComPASS.set_maximum_timestep(250)

communicator = ComPASS.mpi.communicator()
master = ComPASS.mpi.master_proc_rank
rank = ComPASS.mpi.proc_rank

# Collect production temperature on master proc
if rank == master:
    production_temperatures = []
else:
    production_temperatures = None

has_producer = ComPASS.nb_producers() > 0
assert ComPASS.nb_producers() <= 1


def collect_production_temperatures(n, t):
    if has_producer:
        assert ComPASS.nb_producers() == 1
        producers_perforations = list(ComPASS.producers_perforations())
        wellhead_state = producers_perforations[0][
            -1
        ]  # latest perforation corresponds to well head
        production_data = (t, wellhead_state.temperature)
    else:
        production_data = None
    production_data = communicator.gather(production_data, root=master)
    if rank == master:
        production_data = [data for data in production_data if data is not None]
        # We heavily rely on the fact that there is only one production well over the whole field
        # production data can be duplicated has we collect both own and ghost well data
        production_temperatures.append(production_data[0])


#%% First loop: injection of hot water
injection_duration = 0.1 * final_time
standard_loop(
    initial_timestep=1e-5,
    final_time=injection_duration,
    iteration_callbacks=[collect_production_temperatures,],
    output_period=output_period,
)

#%% Second loop: we reinject production water

injectors_data = list(ComPASS.injectors_data())
if injectors_data:
    assert len(injectors_data) == 1
    injector_data = injectors_data[0]
else:
    injector_data = None

# weak coupling: injection at t_n is production at t_{n-1}
def reinject_production(n, t):
    if rank == master:
        production_data = production_temperatures[-1]
    else:
        production_data = None
    production_data = communicator.bcast(production_data, root=master)
    tprod, Tprod = production_data
    assert t == tprod
    if injector_data is not None:
        injector_data.injection_temperature = Tprod


standard_loop(
    initial_timestep=1e-5,
    final_time=final_time,
    iteration_callbacks=[collect_production_temperatures, reinject_production],
    output_period=output_period,
)

if rank == master:
    np.save(
        ComPASS.to_output_directory(
            "production_temperatures_%05d" % ComPASS.mpi.proc_rank
        ),
        production_temperatures,
    )
