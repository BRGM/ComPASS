#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.utils.mesh import extrude

# import vtkwriters as vtkw

p_right = 1.0 * bar  # initial reservoir pressure
T_right = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1e-12  # column permeability in m^2 (low permeability -> bigger time steps)
K_matrix = 2  # bulk thermal conductivity in W/m/K
phi_matrix = 0.15  # column porosity
mass_flux = 1e-1

Lx = 100.0
nx, nz = 10, 3  # number of vertices along axis

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)  # no gravity

mesh = None
if mpi.is_on_master_proc:
    # quads
    assert nx >= 2 and nz >= 2  # 2 points for at least one cell
    x = np.linspace(0, Lx, nx)
    dx = x[1] - x[0]
    assert dx > 0
    pts = []
    for k in range(nz):
        pts.append(np.vstack([x, np.tile(k * dx, nx)]).T)
    pts = np.vstack(pts)
    row = np.vstack([np.array([0, 1, nx + 1, nx]) + k for k in range(nx - 1)])
    quads = np.vstack([row + k * nx for k in range(nz - 1)])
    mesh = extrude(pts, quads, axis=1, thickness=dx)

    # vtkw.write_vtu(vtkw.vtu_doc(mesh.vertices, mesh.cell_nodes), "test")
    # vtkw.write_vtp(vtkw.vtp_doc(mesh.vertices, mesh.face_nodes), "test-faces")


def right_nodes():
    return simulation.global_vertices()[:, 0] >= Lx


simulation.init(
    mesh=mesh,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    set_dirichlet_nodes=right_nodes,
    cell_thermal_conductivity=K_matrix,
)

if mesh is not None:
    del mesh
    mesh = None

for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
]:
    states.context[:] = 2
    states.p[:] = p_right
    states.T[:] = T_right
    states.S[:] = [0, 1]
    states.C[:] = 1.0


Neumann = ComPASS.NeumannBC()
Neumann.molar_flux[:] = mass_flux  # one component
Neumann.heat_flux = mass_flux * simulation.liquid_molar_enthalpy(p_right, T_right)
face_centers = simulation.face_centers()
simulation.set_Neumann_faces(face_centers[:, 0] <= 0, Neumann)

final_time = 1e4 * year
output_period = 0.1 * final_time
standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=TimeStepManager(1 * year, output_period),
    output_period=output_period,
)

simulation.postprocess()
