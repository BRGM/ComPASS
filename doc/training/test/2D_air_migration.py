# Documentation : https://charms.gitlabpages.inria.fr/ComPASS/
import numpy as np
import ComPASS
from ComPASS.utils.units import *  # contains MPa, degC2K, km, year...
from ComPASS.utils.grid import on_zmin, on_zmax
from ComPASS.messages import warning

from fracture_factory import AASFractureNetwork

# -------------------------------------------------------------------
# Geometry
depth = 200
overburden_thickness = 0.2 * depth
nb_fractures = 20
ptop = 1 * bar
T0 = degC2K(20)
gravity = 10.0
seed = 1234  # change this to generate your own mesh !

# -------------------------------------------------------------------
# Reservoir petrophysics
k_matrix = 1e-14  # matrix permeability in m^2
k_overburden = 1e-20  # impermeable overburden
k_fracture = 1e-12  # fracture permeability in m^2
omega_matrix = 0.15  # matrix porosity
omega_fracture = 0.5  # fracture porosity
K = 2  # reservoir thermal conductivity in W/m/K

# -------------------------------------------------------------------
# Set output informations
ComPASS.set_output_directory_and_logfile(__file__)
# Load the water2ph physics : it contains the water component
# which can be in liquid and/or gas phase
simulation = ComPASS.load_physics("immiscible2ph")
simulation.set_gravity(gravity)


# -------------------------------------------------------------------
# Create a Cartesian grid with cubic cells
nx = nz = 20
# nx = nz = 60
ny = 1
dx = dy = dz = depth / nz
Lx = nx * dx
Ly = ny * dy
Lz = nz * dz
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(0.0, 0.0, -depth),
)

# -------------------------------------------------------------------
# Initialize the fracture network setting returning true
# for face centers that are fracture faces
def tag_fracture_faces():
    face_centers = simulation.compute_global_face_centers()
    xc, yc, zc = [face_centers[:, j] for j in range(3)]
    epsilon = min(dx, dy, dz) / 3.0
    assert epsilon > 1e-14
    boundary_faces = (
        (np.abs(xc) < epsilon)
        | (np.abs(xc - Lx) < epsilon)
        | (np.abs(yc) < epsilon)
        | (np.abs(yc - Ly) < epsilon)
        | (np.abs(zc) < epsilon)
        | (np.abs(zc + depth) < epsilon)
    )
    interior_faces = np.logical_not(boundary_faces)
    # no fractures in overburden
    interior_faces[zc >= -overburden_thickness] = False
    # define a fracture network in 2D using x and z oordinates
    dfn = AASFractureNetwork(
        nb_fractures,
        (0, -depth),
        (nx + 1, nz + 1),
        (dx, dz),
        max(nx, nz) // 2,
        seed=seed,
    )
    # keep x, z center coordinates
    points = face_centers[interior_faces][:, (0, 2)]
    fractures = dfn(points, threshold=epsilon)
    interior_faces[interior_faces] = fractures
    if not np.any(interior_faces):
        warning("No fracture face selected !!!")
    return interior_faces


def matrix_permeability():
    # set the overburden permeability
    return something


simulation.init(
    mesh=grid,
    cell_porosity=omega_matrix,
    cell_thermal_conductivity=K,
    cell_permeability=matrix_permeability,
    fracture_faces=tag_fracture_faces,
    fracture_porosity=omega_fracture,
    fracture_thermal_conductivity=K,
    fracture_permeability=k_fracture,
)


# -------------------------------------------------------------------
# Initialize the domain with liquid phase

Xl = simulation.build_state(simulation.Context.diphasic, p=ptop, T=T0, Sg=0)
Xg = simulation.build_state(simulation.Context.diphasic, p=ptop, T=T0, Sg=1)

rho = 1000.0  # water density approximation at 1 bar, 20°C

all_states = simulation.all_states()
all_states.set(Xl)
z = simulation.all_positions()[:, 2]
all_states.p[:] = ptop - rho * gravity * z

# -------------------------------------------------------------------
# Identify and set the Dirichlet nodes

# -------------------------------------------------------------------
# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

lsolver = linear_solver(simulation, direct=True)
# lsolver = linear_solver(simulation, tolerance=1e-10, restart_size=40, max_iterations=400)
newton = Newton(simulation, 1e-5, 8, lsolver)

# -------------------------------------------------------------------
# Execute the time loop with solver parameters
simulation.standard_loop(
    initial_timestep=day,
    final_time=2 * year,
    output_period=10 * day,
    newton=newton,
)


# -------------------------------------------------------------------
# Some postprocesses, it allows to visualize with Paraview
simulation.postprocess()
