import numpy as np

from ComPASS.mpi import MPI


def total_phase_volume(simulation, phase):
    # FIXME: use PartInfo
    nn = simulation.number_of_nodes()
    non = simulation.number_of_own_nodes()
    nf = simulation.number_of_fractures()
    nof = simulation.number_of_own_fractures()
    # nc = simulation.number_of_cells()
    noc = simulation.number_of_own_cells()
    nodes = slice(0, non)
    fractures = slice(nn, nn + nof)
    cells = slice(nn + nf, nn + nf + noc)
    phi = simulation.phase_index(phase)
    states = simulation.all_states()
    Sphi = states.S[:, phi].ravel()
    porous_volume = np.array(simulation.porous_volume_Darcy(), copy=False)
    result = 0
    for locus in [nodes, fractures, cells]:
        result += np.sum(Sphi[locus] * porous_volume[locus])
    return MPI.COMM_WORLD.allreduce(result, MPI.SUM)
