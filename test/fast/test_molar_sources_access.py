import numpy as np
import ComPASS
from ComPASS._kernel import get_kernel


def test_molar_sources_access():

    ComPASS.set_output_directory_and_logfile(__file__)
    simulation = ComPASS.load_physics("diphasic")

    ncomp = 2

    L = 1
    nx, ny, nz = 3, 2, 1

    grid = ComPASS.Grid(
        shape=(nx, ny, nz),
    )

    def u(pts):
        x, y, z = [pts[:, j] for j in range(3)]
        return 0 * np.cos((2 * np.pi / L) * x) * np.sin((2 * np.pi / L) * y)

    def cell_molar_sources():
        centers = simulation.compute_global_cell_centers()
        res = np.zeros((len(centers), ncomp))
        res[:] = 999
        return res
        # res = -((2 * np.pi / L) ** 2) * u(centers)
        # res[:] = 0
        # res[np.linalg.norm(centers, axis=1) < 0.2 * L] = 1
        # return res

    simulation.init(
        mesh=grid,
        cell_porosity=0.1,
        cell_permeability=1.0e-12,
        cell_thermal_conductivity=2.0,
        cell_molar_sources=cell_molar_sources,
    )

    print(simulation.all_molar_sources_vol())
    print(simulation.node_molar_sources_vol())
    print(simulation.cell_molar_sources_vol())
    print(simulation.fracture_molar_sources_vol())
    print(simulation.cell_molar_sources())


if __name__ == "__main__":
    test_molar_sources_access()
