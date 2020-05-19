"""

This module acts like a singleton object holding the state of the simulation.

In a futur evolution, it should become a class that accepts several instances.

"""

import os
import atexit

import numpy as np

import MeshTools as MT
import MeshTools.GridTools as GT
from MeshTools.vtkwriters import vtk_celltype


from .. import mpi
from .. import runtime
from .. import newton
from ..RawMesh import RawMesh
from ..distributed_system import DistributedSystem
from ..ghosts.synchronizer import Synchronizer
from ..petsc import PetscElements
from .. import messages

from .._kernel import get_kernel
from .._kernel import simulation_wrapper as _sw

from .utils import (
    reshape_as_scalar_array, reshape_as_tensor_array
)
from .data import set_fractures, collect_all_edges
from . import state
from . import _simulation_object


# This is temporary but will be generalized in the future
# here Properties will just be used as a namespace
# FIXME never used -> to delete ?
class Properties:
    pass


# FIXME: transient... to be set elsewhere
def default_Newton():
    get_kernel()
    # Legacy parameters
    # NB: you should remove 1 iteration in comparison with legacy values
    assert state.info.petsc is not None
    return newton.Newton(
        _simulation_object.self,
        1e-5, 8, newton.LinearSolver(1e-6, 150),
        solver_fmk=state.info.petsc
    )


# FIXME: grid is kept for backward compatibility, should be deprecated
#        then mesh should not default to None and be explicitely provided
def init(
    grid = None,
    mesh = None,
    wells = lambda: [],
    display_well_ids = False,
    fracture_faces = lambda: None,
    set_dirichlet_nodes = lambda: None,
    set_pressure_dirichlet_nodes = lambda: None,
    set_temperature_dirichlet_nodes = lambda: None,
    set_global_flags = None,
    set_global_rocktype = None,
    mesh_parts = None,
    **kwargs
):
    """
    Initialize many simulation properties and distribute the mesh.
    Before a call to this function, the mesh is global in the sense that
    there is only one mesh on the master proc.
    After its execution the mesh is distributed, i.e.
    there is as many local meshes (with possibly ghost elements) as procs.

    Most of the work will be performed on the master processor.
    At the end of the method the mesh is partitioned and properties are distributed
    to all procs.

    :param mesh: the mesh the simulation is run on, it can be a grid generated with ComPASS.Grid,
        or a MeshTools object. MeshTools provides several helpers modules to load and modify meshes.

    :param wells: a python sequence of well objects
    :param display_well_ids: a boolean to output well ids at the beginning of the simulation (defaults to False)
    :param fracture_faces: the face id of faces that are to be considered as fractures,
        it can also be a mask over all faces
    :param set_dirichlet_nodes: the ids of all nodes that hold boundary conditions (pressure + temperature)
        it can also be a mask over all nodes
    :param set_pressure_dirichlet_nodes: the ids of all nodes that hold constant pressure boundary conditions
        it can also be a mask over all nodes
    :param set_temperature_dirichlet_nodes: the ids of all nodes that hold constant temperature boundary conditions
        it can also be a mask over all nodes

    Petrophysical parameters are mandatory depending on the elements that are present (if there are no fractures,
    fracture properties are not mandatory).

    :param cell_permeability: can be a scalar, the diagonal part or the full permeanility tensor, or an array
        of any of these three types  with as many elements as mesh cells  
    :param cell_porosity: can be a scalar or an array with as many elements as mesh cells
    :param cell_thermal_conductivity: can be a scalar, the diagonal part or the full permeanility tensor, or an array
        of any of these three types with as many elements as mesh cells 
    :param fracture_permeability: can be a scalar or an array with as many elements as fracture cells 
        (fracture permeability is isotropic as doing otherwise would require fracture orientation) 
    :param fracture_porosity: can be a scalar or an array with as many elements as fracture cells
    :param fracture_thermal_conductivity: can be a scalar

    Some parameters control the numerical scheme:

    :param cell_omega_Darcy: the cell volume proportion that is distributed
        to nodes for the discretisation of the pressure gradient (Darcy law)
    :param cell_omega_Fourier: the fracture volume proportion that is distributed
        to nodes for the discretisation of the pressure gradient (Darcy law)
    :param fracture_omega_Darcy: the cell volume proportion that is distributed
        to nodes for the discretisation of the temperature gradient (Fourier law)
    :param fracture_omega_Fourier: the fracture volume proportion that is distributed
        to nodes for the discretisation of the temperature gradient (Fourier law)
    """
    kernel = get_kernel()
    assert not state.mesh_is_local
    # here properties will just be used as a namespace
    properties = {}
    def call_if_callable(f):
        if callable(f):
            return f()
        return f
    def make_property_accessor(value):
        return lambda: value
    for location in ['cell', 'fracture']:
        for property in [
            'porosity', 'permeability', 'thermal_conductivity',
            'omega_Darcy', 'omega_Fourier'
        ]:
            name = '%s_%s' % (location, property)
            try:
                f = kwargs[name]
                if callable(f):
                    properties[name] = f
                else:
                    properties[name] = make_property_accessor(f)
            except KeyError:
                properties[name] = lambda: None
    # FUTURE: This could be managed through a context manager ?
    assert not state.initialized
    # FIXME: grid is kept for backward compatibility, should be deprecated
    if type(mesh) is str:
        if not os.path.exists(mesh):
            print('Mesh file (%s) not found!' % mesh)
        messages.error('Loading mesh from file is desactivated.')
    else:
        if grid is not None:
            messages.deprecation('Use mesh keyword instead of grid')
            if mesh is not None:
                messages.error('You cannot define both grid and mesh keywords')
            mesh = grid
    init_and_load_mesh(mesh)
    if mpi.is_on_master_proc and set_global_flags is not None:
        assert callable(set_global_flags)
        set_global_flags()
    if mpi.is_on_master_proc:
        well_list = list(wells())
        kernel.set_well_geometries(well_list)
        kernel.global_mesh_mesh_bounding_box()
        kernel.global_mesh_compute_all_connectivies()
        check_well_geometry(well_list)
        fractures = call_if_callable(fracture_faces)
        if fractures is not None:
            set_fractures(fractures)
        kernel.global_mesh_node_of_frac()
        _simulation_object.self.set_global_dirichlet_nodes(
            call_if_callable(set_dirichlet_nodes),
            call_if_callable(set_pressure_dirichlet_nodes),
            call_if_callable(set_temperature_dirichlet_nodes)
        )
        kernel.global_mesh_frac_by_node()
        # The following line is necessary to allocate arrays in the fortran code
        kernel.global_mesh_allocate_rocktype()
        if set_global_rocktype is not None:
            assert callable(set_global_rocktype)
            set_global_rocktype()
        _set_petrophysics_on_global_mesh(properties, fractures)
        if 'cell_heat_source' in kwargs:
            value = call_if_callable(kwargs['cell_heat_source'])
            _set_property_on_global_mesh('heat_source', 'cell', value)
        kernel.global_mesh_make_post_read_well_connectivity_and_ip()
        # FIXME: we should distinguish well nature and well status
        for well in well_list:
            if well.closed:
                print("WARNING")
                print("WARNING")
                print("WARNING Closed well will be DISCARDED.")
                print("WARNING Set the well as a producer or injector,")
                print("WARNING and close it before running the simulation.")
                print("WARNING")
                print("WARNING")
                break
        kernel.set_well_data(well_list, display_well_ids)
        kernel.compute_well_indices()
        summarize_simulation()
        if mesh_parts is None:
            mesh_parts = part_mesh()
        else:
            ucolors = np.unique(mesh_parts)
            assert ucolors.min()>=0
            assert ucolors.max()<mpi.communicator().size
        kernel.init_phase2_partition(mesh_parts)
    state.mesh_is_local = True
    mpi.synchronize()
    kernel.init_phase2_build_local_mesh()
    kernel.init_phase2_setup_contexts()
    setup_VAG(properties)
    kernel.init_phase2_setup_solvers(state.info.activate_cpramg, state.info.activate_direct_solver)
    system = DistributedSystem(kernel)
    state.info.system = system
    state.info.ghosts_synchronizer = Synchronizer(system)
    state.info.petsc = PetscElements(system)
    mpi.synchronize() # wait for every process to synchronize
    # FUTURE: This could be managed through a context manager ?
    state.initialized = True
    atexit.register(_exit_eos_and_finalize)


def init_and_load_mesh(mesh):
    kernel = get_kernel()
    kernel.init_warmup(runtime.logfile)
    if mpi.is_on_master_proc:
        # this is a bit redundant but we want to rely entirely on MeshTools
        if type(mesh) is RawMesh:
            subsizes = lambda collection: np.array([len(a) for a in collection])
            make_pointers = lambda a: np.cumsum(np.hstack([[0], a]))
            cell_nbnodes = subsizes(mesh.cell_nodes)
            face_nbnodes = subsizes(mesh.face_nodes)
            int_array = lambda a: np.asarray(a, dtype=np.int32)
            double_array = lambda a: np.asarray(a, dtype=np.double)
            kernel.create_mesh(
                double_array(mesh.vertices),
                int_array(make_pointers(cell_nbnodes)),
                int_array(np.hstack(mesh.cell_nodes)),
                int_array(make_pointers(subsizes(mesh.cell_faces))),
                int_array(np.hstack(mesh.cell_faces)),
                int_array(make_pointers(face_nbnodes)),
                int_array(np.hstack(mesh.face_nodes)),
            )
            # cell and face types default to -1
            # try hint with simple geometries
            celltypes = _sw.global_celltypes()
            celltypes[:] = -1
            celltypes[cell_nbnodes==4] = vtk_celltype['tet']
            celltypes[cell_nbnodes==8] = vtk_celltype['voxel']
            facetypes = _sw.global_facetypes()
            print('>'*10, facetypes.shape)
            facetypes[:] = -1
            facetypes[face_nbnodes==3] = vtk_celltype['triangle']
            facetypes[face_nbnodes==4] = vtk_celltype['pixel']
        else:
            if type(mesh) is GT.GridInfo:
                mesh = MT.grid3D(gridinfo=mesh)
            vertices = MT.as_coordinate_array(mesh.vertices)
            cells_nodes, cells_faces, faces_nodes = mesh.COC_data()
            kernel.create_mesh(vertices,
                                *cells_nodes,
                                *cells_faces,
                                *faces_nodes
                                )
            # distribute cell types to be able to reconstruct local meshes
            celltypes = _sw.global_celltypes()
            celltypes[:] = mesh.cells_vtk_ids().astype(np.int8, copy=False)
            facetypes = _sw.global_facetypes()
            facetypes[:] = mesh.faces_vtk_ids().astype(np.int8, copy=False)


def check_well_geometry(wells):
    if len(wells)==0:
        return
    def to_edge_array(raw_edges):
        assert raw_edges.ndim==2 and raw_edges.shape[1]==2
        edge_type = np.dtype([('v1', np.uint64), ('v2', np.uint64)])
        edges = np.empty(raw_edges.shape[0], dtype=edge_type)
        for i, edge in enumerate(raw_edges):
            edges[i] = tuple(edge)
        return edges
    edges = to_edge_array(collect_all_edges())
    for wi, well in enumerate(wells):
        well_segments = well.geometry.segments
        well_edges = np.unique(np.sort(well_segments, axis=1), axis=0)
        assert well_edges.shape==well_segments.shape, 'duplicate well edge'
        well_edges = to_edge_array(well_edges)
        if not np.all(np.isin(well_edges, edges)):
            message = 'Well %d has edges that are not in mesh edges' % wi
            explanation = 'Well %d edges:\n' % wi + str(well_segments) + '\n'
            messages.error(message, explanation=explanation)


def phase_index(phase):
    assert type(phase) is _sw.Phase
    return int(phase)-1 # Fortran indexing


def context_index(context):
    assert type(context) is _sw.Context
    return int(context)-1 # Fortran indexing


def total_accumulation(reset_states=True):
    kernel = get_kernel()
    # FIXME: reset_states is needed because we also use this function in legacy newton convergence
    if reset_states:
        # Enforce Dirichlet values
        kernel.DirichletContribution_update()
        # Update local jacobian contributions (closure laws)
        kernel.IncPrimSecd_update_secondary_dependencies() # FIXME: this is needed to update globals used in LoisThermoHydro_compute
        kernel.LoisThermoHydro_compute() 
        kernel.Residu_update_accumulation()
    local = np.zeros(_sw.Residuals.npv(), dtype=np.double)
    for states in [
        _sw.own_node_states(),
        _sw.own_fracture_states(),
        _sw.own_cell_states()
    ]:
        local+= np.linalg.norm(states.accumulation, 1, axis=0)
    total = np.zeros(_sw.Residuals.npv(), dtype=np.double)
    mpi.communicator().Allreduce(local, total, mpi.MPI.SUM)
    return total


def _node_info():
    return np.rec.array(_sw.node_info(), copy=False)


def all_fracture_edges_tagged():
    faces = _sw.frac_face_id() - 1 # Fortran indexing...
    face_nodes = _sw.get_connectivity().NodebyFace
    info = _node_info()
    for face in faces:
        nodes = np.array(face_nodes[face], copy=False) - 1 # Fortran indexing...
        in_fracture = info.frac[nodes] == ord('y')
        if not all(in_fracture):
            return False
    return True


def find_fracture_edges(faces):
    if faces.dtype==np.bool:
        faces = np.nonzero(faces)[0]
    face_nodes = _sw.get_connectivity().NodebyFace
    # we do not want to store twice the same edge
    fracture_edges = set()
    info = _node_info()
    for face in faces:
        nodes = np.array(face_nodes[face], copy=False) - 1 # Fortran indexing...
        in_fracture = info.frac[nodes] == ord('y')
        for i in range(nodes.shape[0]):
            if in_fracture[i-1] and in_fracture[i]:
                if nodes[i-1]<nodes[i]:
                    fracture_edges.add((nodes[i-1], nodes[i]))
                else:
                    fracture_edges.add((nodes[i], nodes[i-1]))
    return  np.array(list(fracture_edges))


def pressure_dirichlet_nodes():
    return _node_info().pressure == ord('d')


def pressure_dirichlet_values():
    where = pressure_dirichlet_nodes()
    result = np.tile(np.nan, where.shape)
    result[where] = _sw.dirichlet_node_states().p[where]
    return result


def temperature_dirichlet_nodes():
    return _node_info().temperature == ord('d')


def temperature_dirichlet_values():
    where = temperature_dirichlet_nodes()
    result = np.tile(np.nan, where.shape)
    result[where] = _sw.dirichlet_node_states().T[where]
    return result


def dirichlet_nodes():
    return pressure_dirichlet_nodes() | temperature_dirichlet_nodes()



#---------------------------------------------------------------------------#
# distribution
#


def cell_distribution(colors):
    """Distibute cells into list for each proc according to colors
    and add a laer of ghost cells."""
    assert mpi.is_on_master_proc
    nprocs = mpi.communicator().size
    assert np.all(
        np.in1d(
            np.unique(colors), np.arange(nprocs),
            assume_unique=True
        )
    )
    domains = [
        np.nonzero(colors==color)[0]
        for color in range(nprocs)
    ]
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
        ghost_layers.append({
            color: ghosts[ghosts_color==color]
            for color in np.unique(ghosts_color)
        })
    return domains, ghost_layers


def part_mesh():
    kernel = get_kernel()
    nparts = mpi.communicator().size
    if nparts>1:
        cell_colors = kernel.metis_part_graph(
            _sw.get_global_connectivity().CellbyCell, nparts,
        )
        cell_colors-= 1 # Fortran indexing
    else:
        cell_colors = np.zeros(_sw.global_number_of_cells(), dtype=np.int32)
    return cell_colors


#---------------------------------------------------------------------------#
# simulation output
#


@mpi.on_master_proc
def summarize_simulation():
    kernel = get_kernel()
    print('  NbCell:      ', _sw.global_number_of_cells())
    # print('  NbFace:      ', NbFace)
    print('  NbNode:      ', _sw.global_number_of_nodes())
    # print('  NbFrac:      ', NbFrac)
    # print('  NbWellInj    ', NbWellInj)
    # print('  NbWellProd   ', NbWellProd)
    # print('  NbDirNode P: ', NbDirNodeP)
    if kernel.has_energy_transfer_enabled:
        print('Energy transfer enabled')
        # print('  NbDirNode T: ', NbDirNodeT)
    print('  Ncpus :     ', mpi.communicator().size)


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




#---------------------------------------------------------------------------#
# Internal
#


def setup_VAG(properties):
    kernel = get_kernel()
    omegas = {}
    for law in ('Darcy', 'Fourier'):
        for location in ('cell', 'fracture'):
            omega = properties[location + '_omega_' + law]()
            if omega is None:
                omega = 0.075 if location=='cell' else 0.15
            omegas['%s_%s'%(location, law)] = omega
    kernel.VAGFrac_TransDarcy()
    if kernel.has_energy_transfer_enabled():
        kernel.VAGFrac_TransFourier()
    kernel.VAGFrac_VolsDarcy(omegas['cell_Darcy'], omegas['fracture_Darcy'])
    if kernel.has_energy_transfer_enabled():
        kernel.VAGFrac_VolsFourier(omegas['cell_Fourier'], omegas['fracture_Fourier'])


def _set_property_on_global_mesh(property, location, value, fractures=None):
    kernel = get_kernel()
    buffer = np.array(
        getattr(kernel, 'get_global_%s_%s_buffer' % (location, property))(),
        copy = False,
    )
    n = buffer.shape[0]
    dim = 3
    if location=='fracture':
        assert fractures is not None, (
            "no fractures while setting %s on %s" %(property, location)
        )            
        n = np.count_nonzero(fractures)
        dim = 2
    if property in ['permeability', 'thermal_conductivity'] and location=='cell':
        value = reshape_as_tensor_array(value, n, dim)
    else:
        value = reshape_as_scalar_array(value, n)
    if location=='fracture':
        buffer[fractures] = value
    else:
        buffer[:] = value


def _set_petrophysics_on_global_mesh(properties, fractures):
    kernel = get_kernel()
    if mpi.is_on_master_proc:
        kernel.global_mesh_allocate_petrophysics()
        useful_properties =  ['porosity', 'permeability']
        if kernel.has_energy_transfer_enabled():
            useful_properties.append('thermal_conductivity')
        for location in ['cell', 'fracture']:
            for property in useful_properties:
                value = properties[location + '_' + property]()
                if value is not None:
                    _set_property_on_global_mesh(
                        property, location, value, fractures
                    )
                else:
                    if location=='cell':
                        messages.error('You must define: cell_%s' % property)
                    elif fractures is not None:
                        messages.error('You must define: fracture_%s' % property)
        print('petrophysics')
        _petrophysics_statistics_on_global_mesh(fractures)
        kernel.global_mesh_set_all_rocktypes()


def _petrophysics_statistics_on_global_mesh(fractures):
    kernel = get_kernel()
    if mpi.is_on_master_proc:
        for location in ['cell', 'fracture']:
            # TODO permeability and thermal_condutvity are tensors
            # for property in ['porosity', 'permeability', 'thermal_conductivity']:
            for property in ['porosity']:
                buffer = np.array(getattr(kernel, 'get_global_%s_%s_buffer' % (location, property))(), copy = False)
                if location=='fracture':
                    buffer = buffer[fractures]
                print(buffer.min(), '<=', location, property, '<=', buffer.max())


def _exit_eos_and_finalize():
    _sw.finalize_model()
    _sw.finalize()
