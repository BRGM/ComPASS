"""

This module defines a class that is a singleton object holding the state of the simulation.

In a futur evolution, the singleton constraint should be released.

"""

import os
import atexit
import numpy as np

from .. import mpi
from ..distributed_system import DistributedSystem
from ..ghosts.synchronizer import Synchronizer
from .. import messages

from .data import set_fractures
from .base import (
    init_and_load_mesh,
    _set_petrophysics_on_global_mesh,
    _set_property_on_global_mesh,
    summarize_simulation,
    dump_global_mesh,
    part_mesh,
    setup_scheme,
    _exit_physics_and_finalize,
)

# FIXME: grid is kept for backward compatibility, should be deprecated
#        then mesh should not default to None and be explicitely provided
def init(
    simulation,
    grid=None,
    mesh=None,
    cell_porosity=lambda: None,
    cell_permeability=lambda: None,
    cell_thermal_conductivity=lambda: None,
    cell_omega_Darcy=lambda: None,
    cell_omega_Fourier=lambda: None,
    fracture_porosity=lambda: None,
    fracture_permeability=lambda: None,
    fracture_thermal_conductivity=lambda: None,
    fracture_omega_Darcy=lambda: None,
    fracture_omega_Fourier=lambda: None,
    cell_heat_source=None,
    wells=lambda: [],
    display_well_ids=False,
    fracture_faces=lambda: None,
    set_dirichlet_nodes=lambda: None,
    set_pressure_dirichlet_nodes=lambda: None,
    set_temperature_dirichlet_nodes=lambda: None,
    set_global_flags=None,
    set_global_rocktype=None,
    mesh_parts=None,
    use_Kway_part_graph=False,
    well_model="single_phase",
    dump_mesh_before_distribution=False,
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
    :param mesh_parts: directly provide mesh partitioning with a sequence of integer representing a color
    :param use_Kway_part_graph: if True use METIS_PartGraphKway to part mesh else use METIS_PartGraphRecursive (the default)

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
    :param cell_heat_source: volumic thermal source term (will be multiplicated by the cell volume)
        for each cell

    Petrophysical parameters are mandatory depending on the elements that are present (if there are no fractures,
    fracture properties are not mandatory).

    :param cell_permeability: can be a scalar, the diagonal part or the full permeability tensor, or an array
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
    :param well_model: the well model to be used (can be "single_phase" (default) or "two_phases")

    You can dump the global mesh before partitioning and properties distribution using
    the `dump_mesh_before_distribution` option which defaults to `False`.

    :param dump_mesh_before_distribution: boolean value that defaults to `False`

    """
    kernel = simulation.get_kernel()
    assert not simulation.mesh_is_local
    # here properties will just be used as a namespace
    properties = {}

    def call_if_callable(f):
        if callable(f):
            return f()
        return f

    def make_property_accessor(value):
        return lambda: value

    for location in ["cell", "fracture"]:
        for property in [
            "porosity",
            "permeability",
            "thermal_conductivity",
            "omega_Darcy",
            "omega_Fourier",
        ]:
            name = "%s_%s" % (location, property)
            f = eval(name)
            if callable(f):
                properties[name] = f
            else:
                properties[name] = make_property_accessor(f)
    # FUTURE: This could be managed through a context manager ?
    assert not simulation.initialized
    # FIXME: grid is kept for backward compatibility, should be deprecated
    if type(mesh) is str:
        if not os.path.exists(mesh):
            messages.warning("Mesh file (%s) not found!" % mesh)
        messages.error("Loading mesh from file is desactivated.")
    else:
        if grid is not None:
            messages.deprecation("Use mesh keyword instead of grid")
            if mesh is not None:
                messages.error("You cannot define both grid and mesh keywords")
            mesh = grid
    init_and_load_mesh(mesh)
    if mpi.is_on_master_proc:
        well_list = list(wells())
        kernel.set_well_geometries(well_list)
        kernel.global_mesh_mesh_bounding_box()
        kernel.global_mesh_compute_all_connectivies()
        if set_global_flags is not None:
            assert callable(set_global_flags)
            set_global_flags()
        simulation.check_well_geometry(well_list)
        fractures = call_if_callable(fracture_faces)
        if fractures is not None:
            set_fractures(fractures)
        kernel.global_mesh_node_of_frac()
        simulation.set_global_dirichlet_nodes(
            call_if_callable(set_dirichlet_nodes),
            call_if_callable(set_pressure_dirichlet_nodes),
            call_if_callable(set_temperature_dirichlet_nodes),
        )
        kernel.global_mesh_frac_by_node()
        # The following line is necessary to allocate arrays in the fortran code
        kernel.global_mesh_allocate_rocktype()
        if set_global_rocktype is not None:
            assert callable(set_global_rocktype)
            set_global_rocktype()
        _set_petrophysics_on_global_mesh(properties, fractures)
        if cell_heat_source is not None:
            value = call_if_callable(cell_heat_source)
            _set_property_on_global_mesh("heat_source", "cell", value)
        kernel.build_well_connectivity()
        # FIXME: we should distinguish well nature and well status
        for well in well_list:
            if well.closed:
                messages.warning(
                    "Closed well will be DISCARDED.\n"
                    "Set the well as a producer or injector "
                    "and close it before running the simulation.",
                )
                break
        kernel.set_well_data(well_list, display_well_ids)
        kernel.compute_well_indices()
        summarize_simulation()
        if dump_mesh_before_distribution:
            dump_global_mesh(simulation)
        part_mesh(simulation, mesh_parts=mesh_parts, use_Kway=use_Kway_part_graph)
    simulation.mesh_is_local = True
    mpi.synchronize()
    kernel.init_phase2_build_local_mesh()
    kernel.init_phase2_setup_contexts()
    setup_scheme(simulation, properties)
    kernel.init_phase2_setup_solvers()
    system = DistributedSystem(kernel)
    simulation.info.system = system
    simulation.info.ghosts_synchronizer = Synchronizer(system)
    if kernel.has_freeflow_structures:  # FIXME: to be removed
        kernel.clear_freeflow_faces()
    mpi.synchronize()  # wait for every process to synchronize
    if kernel.get_fill_kr_arrays() is not None:
        messages.warning("not overriding kr function in init with default value")
    else:
        simulation.set_kr_functions()
    if kernel.get_fill_phase_pressure_arrays() is not None:
        messages.warning("not overriding Pc function in init with default value")
    else:
        simulation.set_phase_pressure_functions()
    simulation.set_well_model(well_model)
    assert simulation.unknown_producers_density
    # FUTURE: This could be managed through a context manager ?
    simulation.initialized = True
    atexit.register(_exit_physics_and_finalize)
