"""

This module acts like a singleton object holding the state of the simulation.

In a futur evolution, it should become a class that accepts several instances.

"""

import os
from pathlib import Path

import numpy as np

import MeshTools as MT
import MeshTools.GridTools as GT
import vtkwriters as vtkw

from ComPASS import metis
from .. import mpi
from .. import runtime
from ..RawMesh import RawMesh
from .. import messages
from ..utils.adjacencies import filter_adjacency_table

from .._kernel import get_kernel
from .._kernel import simulation_wrapper as _sw

from ..schemes.VAG import VAGScheme
from . import utils

# This is temporary but will be generalized in the future
# here Properties will just be used as a namespace
# FIXME never used -> to delete ?
class Properties:
    pass


def init_and_load_mesh(mesh):
    kernel = get_kernel()
    kernel.init_warmup(runtime.logfile)
    if mpi.is_on_master_proc:
        # this is a bit redundant but we want to rely entirely on MeshTools
        if isinstance(mesh, RawMesh):
            cn = mesh.get_cell_nodes()
            cf = mesh.get_cell_faces()
            fn = mesh.get_face_nodes()

            kernel.create_mesh(
                mesh.get_vertices(),
                cn.pointers,
                cn.contents,
                cf.pointers,
                cf.contents,
                fn.pointers,
                fn.contents,
            )

            mesh.fill_cell_types(_sw.global_celltypes())
            mesh.fill_face_types(_sw.global_facetypes())
        else:
            if type(mesh) is GT.GridInfo:
                mesh = MT.grid3D(gridinfo=mesh)
            vertices = MT.as_coordinate_array(mesh.vertices)
            cells_nodes, cells_faces, faces_nodes = mesh.COC_data()
            kernel.create_mesh(vertices, *cells_nodes, *cells_faces, *faces_nodes)
            # distribute cell types to be able to reconstruct local meshes
            celltypes = _sw.global_celltypes()
            celltypes[:] = mesh.cells_vtk_ids().astype(np.int8, copy=False)
            facetypes = _sw.global_facetypes()
            facetypes[:] = mesh.faces_vtk_ids().astype(np.int8, copy=False)


def phase_index(phase):
    assert type(phase) is _sw.Phase
    return int(phase) - 1  # Fortran indexing


def context_index(context):
    assert type(context) is _sw.Context
    return int(context) - 1  # Fortran indexing


def total_accumulation(reset_states=True):
    kernel = get_kernel()
    # FIXME: reset_states is needed because we also use this function in legacy newton convergence
    if reset_states:
        kernel.LoisThermoHydro_compute_phase_pressures()  # phase pressures are needed in IncPrimSecd_update_secondary_dependencies (fugacities)
        # Enforce Dirichlet values
        kernel.DirichletContribution_update()
        # Update local jacobian contributions (closure laws)
        kernel.IncPrimSecd_update_secondary_dependencies()  # FIXME: this is needed to update globals used in LoisThermoHydro_compute
        kernel.LoisThermoHydro_compute_phase_pressures_derivatives()  # needs to be done IncPrimSecd_update_secondary_dependencies (needs NumIncTotalPrim)
        kernel.LoisThermoHydro_compute()
        kernel.Residu_update_accumulation()
    local = np.zeros(_sw.Residuals.npv(), dtype=np.double)
    for states in [
        _sw.own_node_states(),
        _sw.own_fracture_states(),
        _sw.own_cell_states(),
    ]:
        local += np.linalg.norm(states.accumulation, 1, axis=0)
    total = np.zeros(_sw.Residuals.npv(), dtype=np.double)
    mpi.communicator().Allreduce(local, total, mpi.MPI.SUM)
    return total


def _node_info():
    return np.rec.array(_sw.node_info(), copy=False)


def all_fracture_edges_tagged():
    faces = _sw.frac_face_id() - 1  # Fortran indexing...
    face_nodes = _sw.get_connectivity().NodebyFace
    info = _node_info()
    for face in faces:
        nodes = np.array(face_nodes[face], copy=False) - 1  # Fortran indexing...
        in_fracture = info.frac[nodes] == ord("y")
        if not all(in_fracture):
            return False
    return True


def find_fracture_edges(faces):
    if faces.dtype == bool:
        faces = np.nonzero(faces)[0]
    face_nodes = _sw.get_connectivity().NodebyFace
    # we do not want to store twice the same edge
    fracture_edges = set()
    info = _node_info()
    for face in faces:
        nodes = np.array(face_nodes[face], copy=False) - 1  # Fortran indexing...
        in_fracture = info.frac[nodes] == ord("y")
        for i in range(nodes.shape[0]):
            if in_fracture[i - 1] and in_fracture[i]:
                if nodes[i - 1] < nodes[i]:
                    fracture_edges.add((nodes[i - 1], nodes[i]))
                else:
                    fracture_edges.add((nodes[i], nodes[i - 1]))
    return np.array(list(fracture_edges))


def pressure_dirichlet_nodes():
    return _node_info().pressure == ord("d")


def pressure_dirichlet_values():
    where = pressure_dirichlet_nodes()
    result = np.tile(np.nan, where.shape)
    result[where] = _sw.dirichlet_node_states().p[where]
    return result


def temperature_dirichlet_nodes():
    return _node_info().temperature == ord("d")


def temperature_dirichlet_values():
    where = temperature_dirichlet_nodes()
    result = np.tile(np.nan, where.shape)
    result[where] = _sw.dirichlet_node_states().T[where]
    return result


def dirichlet_nodes():
    return pressure_dirichlet_nodes() | temperature_dirichlet_nodes()


# ---------------------------------------------------------------------------#
# distribution
#


def cell_distribution(colors):
    """Distibute cells into list for each proc according to colors
    and add a laer of ghost cells."""
    assert mpi.is_on_master_proc
    nprocs = mpi.communicator().size
    assert np.all(np.in1d(np.unique(colors), np.arange(nprocs), assume_unique=True))
    domains = [np.nonzero(colors == color)[0] for color in range(nprocs)]
    # add one layer of ghost cells
    ghost_layers = []
    neighbors = _sw.get_global_connectivity().CellbyCell
    for color, domain in enumerate(domains):
        # FIXME: Fortran indexing in connectivity -> -1
        extended_domain = np.hstack([neighbors[k] for k in domain]) - 1
        extended_domain = np.unique(np.hstack([domain, extended_domain]))
        own = np.in1d(extended_domain, domain, assume_unique=True)
        ghosts = extended_domain[np.logical_not(own)]
        ghosts_color = colors[ghosts]
        ghost_layers.append(
            {color: ghosts[ghosts_color == color] for color in np.unique(ghosts_color)}
        )
    return domains, ghost_layers


def _part_mesh(use_Kway, connectivity_file=None):
    assert mpi.is_on_master_proc
    kernel = get_kernel()
    nparts = mpi.communicator().size
    cell_connectivity = _sw.get_global_connectivity().CellbyCell
    offsets = cell_connectivity.offsets()
    # FIXME: the following if clause is linked to a bug in single-cell mesh connectivity
    if len(offsets) == 2 and np.all(offsets == 0):
        assert _sw.global_number_of_cells() == 1
        return np.zeros(_sw.global_number_of_cells(), dtype=np.int32)
    assert not np.all(offsets == 0)
    neighbors = cell_connectivity.contiguous_content()
    if connectivity_file is not None:
        connectivity_file = Path(connectivity_file)
        connectivity_file.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            connectivity_file,
            neighbors=neighbors,
            offsets=offsets,
        )
    if nparts > 1:
        # Fortran indexing of neighbors
        cell_colors = metis.part_graph(neighbors - 1, offsets, nparts, Kway=use_Kway)
    else:
        cell_colors = np.zeros(_sw.global_number_of_cells(), dtype=np.int32)
    return cell_colors


# def part_mesh(simulation, use_Kway, connectivity_file=None):
def part_mesh(simulation, mesh_parts=None, use_Kway=False):
    assert mpi.is_on_master_proc, "Mesh partioning is assumed to run on master proc."
    if simulation.is_sequential or simulation.global_number_of_cells() == 1:
        if mesh_parts is None:
            mesh_parts = np.tile(
                mpi.master_proc_rank, simulation.global_number_of_cells()
            )
        else:
            assert tuple(np.unique(mesh_parts)) == (mpi.master_proc_rank,)
            pass
    else:
        if mesh_parts is None:
            mesh_parts = _part_mesh(
                use_Kway=use_Kway,
                # connectivity_file=simulation.runtime.to_output_directory(
                #     "mesh/connectivity"
                # ),
            )
    ucolors = np.unique(mesh_parts)
    assert ucolors.min() >= 0
    assert ucolors.max() < mpi.communicator().size
    parts_file = Path(simulation.runtime.to_output_directory("mesh/parts"))
    parts_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(parts_file, mesh_parts)
    get_kernel().init_phase2_partition(mesh_parts)


# ---------------------------------------------------------------------------#
# simulation output
#


@mpi.on_master_proc
def summarize_simulation():
    kernel = get_kernel()
    print("  NbCell:      ", _sw.global_number_of_cells())
    # print('  NbFace:      ', NbFace)
    print("  NbNode:      ", _sw.global_number_of_nodes())
    # print('  NbFrac:      ', NbFrac)
    # print('  NbWellInj    ', NbWellInj)
    # print('  NbWellProd   ', NbWellProd)
    # print('  NbDirNode P: ', NbDirNodeP)
    if kernel.has_energy_transfer_enabled:
        print("Energy transfer enabled")
        # print('  NbDirNode T: ', NbDirNodeT)
    print("  Ncpus :     ", mpi.communicator().size)


def dump_global_mesh(simulation):
    assert mpi.is_on_master_proc
    output_directory = Path(runtime.output_directory) / "before_distribution"
    output_directory.mkdir(parents=True, exist_ok=True)
    simulation.write_polyhedra_vtu_mesh(str(output_directory / "mesh"))
    connectivity = simulation.get_global_connectivity()
    boundaries = utils.get_boundary_faces(connectivity)
    boundaries = utils.get_faces_nodes(connectivity, boundaries)
    boundary_nodes, boundary_faces = filter_adjacency_table(boundaries)
    vertices = simulation.global_vertices()
    vtkw.write_vtp(
        vtkw.vtp_doc(
            vertices[boundary_nodes],
            boundary_faces,
        ),
        str(output_directory / "boundary_faces"),
    )


def output_visualization_files(iteration):
    kernel = get_kernel()
    kernel.output_visu(iteration, runtime.output_directory)


def add_output_subdirectory(path):
    target = runtime.to_output_directory(path)
    if mpi.is_on_master_proc:
        if not os.path.exists(target):
            os.makedirs(target)
    mpi.synchronize()
    assert os.path.isdir(target)


# ---------------------------------------------------------------------------#
# Internal
#


def setup_VAG(properties):
    def get_omega_values(law):
        omega = {}
        for location in ("cell", "fracture"):
            value = properties[location + "_omega_" + law]()
            if value is not None:
                omega[location] = value
        return omega

    return VAGScheme(
        **{f"omega_{law}": get_omega_values(law) for law in ("Darcy", "Fourier")}
    )


def setup_scheme(simulation, properties):
    assert simulation.scheme is None
    simulation.scheme = vag = setup_VAG(properties)
    vag.compute_transmissivities()
    vag.compute_volumes()


# FIXME: we need fractures as argument because global_fracture_property_buffer
#        has the size of faces and we need to pick only fracture faces
def _set_property_on_global_mesh(property, location, value, fractures=None):
    assert location == "cell" or location == "fracture"
    kernel = get_kernel()
    buffer = np.array(
        getattr(kernel, "get_global_%s_%s_buffer" % (location, property))(),
        copy=False,
    )
    n = buffer.shape[0]
    dim = 3
    if location == "fracture":
        assert fractures is not None, "no fractures while setting %s on %s" % (
            property,
            location,
        )
        if fractures.dtype == bool:
            n = np.count_nonzero(fractures)
        else:
            assert fractures.ndim == 1
            n = fractures.shape[0]
        dim = 2

    if property in ["permeability", "thermal_conductivity"] and location == "cell":
        value = utils.reshape_as_tensor_array(value, n, dim)
    elif property == "molar_sources":
        assert (
            value.shape == buffer.shape
        ), f"molar_sources shape must be (NbCell, NbComp), not {value.shape}"
    else:
        value = utils.reshape_as_scalar_array(value, n)
    if location == "fracture":
        buffer[fractures] = value
    else:
        buffer[:] = value


def _set_petrophysics_on_global_mesh(properties, fractures):
    kernel = get_kernel()
    if mpi.is_on_master_proc:
        kernel.global_mesh_allocate_petrophysics()
        useful_properties = ["porosity", "permeability"]
        if kernel.has_energy_transfer_enabled():
            useful_properties.append("thermal_conductivity")
        for location in ["cell", "fracture"]:
            for property in useful_properties:
                value = properties[location + "_" + property]()
                if value is not None:
                    _set_property_on_global_mesh(property, location, value, fractures)
                else:
                    if location == "cell":
                        messages.error("You must define: cell_%s" % property)
                    elif fractures is not None:
                        messages.error("You must define: fracture_%s" % property)
        print("petrophysics")
        _petrophysics_statistics_on_global_mesh(fractures)
        kernel.global_mesh_set_all_rocktypes()


def _petrophysics_statistics_on_global_mesh(fractures):
    kernel = get_kernel()
    if mpi.is_on_master_proc:
        for location in ["cell", "fracture"]:
            # TODO permeability and thermal_condutvity are tensors
            # for property in ['porosity', 'permeability', 'thermal_conductivity']:
            for property in ["porosity"]:
                buffer = np.array(
                    getattr(kernel, "get_global_%s_%s_buffer" % (location, property))(),
                    copy=False,
                )
                if location == "fracture":
                    buffer = buffer[fractures]
                print(buffer.min(), "<=", location, property, "<=", buffer.max())


def _exit_physics_and_finalize():
    get_kernel().release_physics()
    get_kernel().release_simulation()
    _sw.finalize_model()
    _sw.finalize()
