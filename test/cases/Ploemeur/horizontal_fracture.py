import numpy as np

import MeshTools as MT

import ComPASS
from ComPASS.timeloops import standard_loop
from ComPASS.utils.units import *
from ComPASS.utils.wells import create_vertical_well

#%% mesh

def geometric(x0, xn, a, n):
    assert a>=1 and n>0
    dx = np.cumprod(np.hstack([1, np.tile(a, max(n-1, 0))]))
    assert np.all(dx>0)
    x = np.cumsum(dx)
    assert np.all(x[1:]>=x[:-1])
    return x0 + ((xn - x0)/(x[-1]-x[0])) * (x - x[0])

# grid is centered on the middle of the segment between injection and production points
# hence wells positions are (-0.5 * interwell_distance, 0, 0) and (0.5 * interwell_distance, 0, 0)
interwell_distance = 5 # meters
grid_length = 3 * interwell_distance # meters
grid_depth = 2 * interwell_distance # meters
grid_heigth = 1.5 * interwell_distance # meters
assert 0 < interwell_distance
assert interwell_distance < grid_length
nb_steps_x_outside_wells = 5 
nb_steps_y = 5 
nb_steps_z = 5 
steps_x_ratio_outside_wells = 1.3 
steps_y_ratio = 1.3 
steps_z_ratio = 1.5 

half_interwell_distance = 0.5 * interwell_distance
xout = geometric(half_interwell_distance, 0.5 * grid_length,
                     steps_x_ratio_outside_wells, nb_steps_x_outside_wells)
dx = xout[1] - xout[0]
xbetween = np.linspace(-half_interwell_distance, half_interwell_distance, np.ceil(interwell_distance/dx)+1)
x = np.hstack([
    (-xout[1:])[::-1],
    xbetween,
    xout[1:]
])
ypos = geometric(0, 0.5 * grid_depth, steps_y_ratio, nb_steps_y)
y = np.hstack([(-ypos[1:])[::-1], ypos])
zpos = geometric(0, 0.5 * grid_heigth, steps_z_ratio, nb_steps_z)
dz = zpos[1] - zpos[0]
z = np.hstack([(-zpos[1:])[::-1], zpos])

mesh = MT.grid3D(steps=(x, y, z))

MT.to_vtu(mesh, 'ploemeur_mesh.vtu')

#%% simulation

rhof = 1E3               # specific mass in kg/m^3
cpf = 4200               # specific heat in J/kg/K
rhofcpf = rhof * cpf     # volumetric heat capacity
muf = 1E-3               # viscosity Pa.s
pres = 0. * MPa                  # initial reservoir pressure
Tres = degC2K( 10. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 15. )        # injection temperature - convert Celsius to Kelvin degrees
Qm = 1. * ton / hour            # production/injection flowrate
k_fracture = 1E-12               # fracture permeability in m^2
omega_fracture = 0.55            # fracture porosity
K_fracture = 2                                # bulk thermal conductivity in W/m/K
k_matrix = 1E-20            # matrix permeability in m^2
omega_matrix = 0.15            # matrix porosity
K_matrix = 2                                # bulk thermal conductivity in W/m/K

ComPASS.load_eos('linear_water')
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.volumetric_heat_capacity = rhofcpf
fluid_properties.dynamic_viscosity = muf

ComPASS.set_gravity(0)
ComPASS.set_output_directory_and_logfile(__file__)

def build_wells():
    producer = create_vertical_well((half_interwell_distance, 0), well_radius = 0.1, zmin=-dz, zmax=dz)
    producer.operate_on_flowrate = Qm, -1E99 # mass flow rate, pressure limit (minimum)
    producer.produce()
    injector = create_vertical_well((-half_interwell_distance, 0), well_radius = 0.1, zmin=-dz, zmax=dz)
    injector.operate_on_flowrate = -Qm * ton / hour, 1E99 # mass flow rate, pressure limit (maximum)
    injector.inject(degC2K(Tinjection))
    return (producer, injector)

def select_dirichlet_nodes():
    vertices = ComPASS.global_vertices()
    return (vertices[:, 2] == z.min()) | (vertices[:, 2] == z.max()) 

def select_fractures():
    face_centers = ComPASS.compute_global_face_centers()
    return face_centers[:, 2] == 0

ComPASS.init(
    mesh = mesh,
    wells = build_wells,
    cell_porosity = omega_matrix,
    cell_permeability = k_matrix,
    cell_thermal_conductivity = K_matrix,
    fracture_faces = select_fractures,
    fracture_porosity = omega_fracture,
    fracture_permeability = k_fracture,
    fracture_thermal_conductivity = K_fracture,
    set_dirichlet_nodes = select_dirichlet_nodes
)

def set_states(states):
    states.context[:] = 1
    states.p[:] = pres
    states.T[:] = Tres
    states.S[:] = 1.
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
                ComPASS.node_states(),
                ComPASS.fracture_states(),
                ComPASS.cell_states()]:
    set_states(states)

final_time = 2E4
output_period = 0.05 * final_time
ComPASS.set_maximum_timestep( 0.5 * output_period)
standard_loop(initial_timestep = 1E-5, final_time = final_time, output_period = output_period)
