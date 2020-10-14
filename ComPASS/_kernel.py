import importlib

from . import mpi

kernel = None


def get_kernel():
    assert kernel is not None, "kernel not loaded"
    return kernel


def load_eos(eosname):
    global kernel
    assert kernel is None
    kernel = importlib.import_module("ComPASS.eos.%s" % eosname)
    # CHECKME: we replace the behavior: from import ComPASS.eos.eosname import *
    #          there might be a more elegant way to do this
    # gdict = globals()
    # kdict = vars(kernel)
    # for key in kdict:
    #     if key not in gdict:
    #         gdict[key] = kdict[key]
    kernel.init_model()
    from .simulation import _simulation_object

    kernel.register_simulation(_simulation_object.self)
    return _simulation_object.self


class Wrapper:
    def __init__(self, exposed):
        self._exposed = tuple(exposed)

    def __getattr__(self, name):
        if name in self._exposed:
            return getattr(get_kernel(), name)
        else:
            raise AttributeError(f"wrapper has no attribute {name!r}")

    def __dir__(self):
        res = super().__dir__()
        res.extend(self._exposed)
        return res


common_wrapper = Wrapper(
    [
        # types
        "COCcontainer",
        "COCiterator",
        "COC",
        "States",
        "NeumannContributions",
        "NeumannBC",
        "MeshConnectivity",
        "Residuals",
        "PartElement",
        "PartInfo",
        "NewtonIncrements",
        "CTVector",
        "WellGeometry",
        "Well",
        "WellData",
        "PerforationData",
        "PerforationState",
        "WellGroup",
        "WellInformation",
        "WellPerforations",
        "FluidProperties",
        "BlockMatrix",
    ]
)


simulation_wrapper = Wrapper(
    [
        "LinearSystemBuilder",
        "init_model",
        "finalize_model",
        "finalize",
        "set_well_geometries",
        "set_well_data",
        "get_gravity",
        "set_gravity",
        "global_celltypes",
        "global_facetypes",
        "global_node_info",
        "has_energy_transfer_enabled",
        "global_number_of_cells",
        "global_number_of_nodes",
        "dirichlet_node_states",
        "node_states",
        "cell_states",
        "fracture_states",
        "Residuals",
        "production_whp",
        "injection_whp",
        "nb_cells_own",
        "nb_nodes_own",
        "nb_faces_own",
        "nb_fractures_own",
        "get_connectivity",
        "get_nodes_by_fractures",
        "frac_face_id",
        "facetypes",
        "vertices",
        "celltypes",
        "number_of_phases",
        "State",
        "Phase",
        "Component",
        "Context",
        "number_of_components",
        "mass_fluxes",
        "create_mesh",
        "debug_utils_dump_mesh_info",
        "global_mesh_make_post_read",
        "global_mesh_make_post_read_fracture_and_dirBC",
        "global_mesh_allocate_petrophysics",
        "global_mesh_set_all_rocktypes",
        "global_mesh_make_post_read_well_connectivity_and_ip",
        "global_mesh_mesh_bounding_box",
        "global_mesh_compute_all_connectivies",
        "global_mesh_set_frac",
        "global_mesh_node_of_frac",
        "global_mesh_set_dir_BC",
        "global_mesh_frac_by_node",
        "global_mesh_allocate_rocktype",
        "global_mesh_count_dirichlet_nodes",
        "compute_well_indices",
        "set_peaceman_index_threshold",
        "unset_peaceman_index_threshold",
        "get_atm_pressure",
        "set_atm_pressure",
        "get_atm_temperature",
        "set_atm_temperature",
        "get_atm_flux_radiation",
        "set_atm_flux_radiation",
        "get_soil_emissivity",
        "set_soil_emissivity",
        "get_atm_rain_flux",
        "set_atm_rain_flux",
        "get_rock_volumetric_heat_capacity",
        "set_rock_volumetric_heat_capacity",
        "get_fracture_thickness",
        "set_fracture_thickness",
        "all_states",
        "own_dirichlet_node_states",
        "own_node_states",
        "own_fracture_states",
        "own_cell_states",
        "set_Neumann_faces",
        "set_Neumann_fracture_edges",
        "get_global_connectivity",
        "number_of_nodes",
        "number_of_own_nodes",
        "number_of_cells",
        "number_of_own_cells",
        "number_of_faces",
        "number_of_own_faces",
        "number_of_fractures",
        "number_of_own_fractures",
        "number_of_own_injectors",
        "number_of_own_producers",
        "injection_whp",
        "production_whp",
        "global_nodeflags",
        "global_cellflags",
        "global_faceflags",
        "global_celltypes",
        "global_facetypes",
        "nodeflags",
        "cellflags",
        "faceflags",
        "celltypes",
        "facetypes",
        "face_frac_id",
        "frac_face_id",
        "nb_cells_own",
        "nb_faces_own",
        "nb_nodes_own",
        "nb_fractures_own",
        "nb_wellinj_own",
        "nb_wellprod_own",
        "all_thermal_sources",
        "cellthermalsource",
        "nodethermalsource",
        "fracturethermalsource",
        "all_Fourier_porous_volumes",
        "porovolfouriercell",
        "porovolfouriernode",
        "global_node_info",
        "node_info",
        "global_vertices",
        "vertices",
        "cell_centers",
        "face_centers",
        "global_cell_rocktypes",
        "global_node_rocktypes",
        "global_fracture_rocktypes",
        "cell_rocktypes",
        "node_rocktypes",
        "fracture_rocktypes",
        "get_global_id_faces_buffer",
        "get_cell_heat_source_buffer",
        "get_global_cell_porosity_buffer",
        "get_global_fracture_porosity_buffer",
        "get_global_cell_permeability_buffer",
        "get_global_fracture_permeability_buffer",
        "get_global_cell_thermal_conductivity_buffer",
        "get_global_fracture_thermal_conductivity_buffer",
        "porous_volume_Darcy",
        "petrophysics",
        "set_fill_kr_arrays",
        # "cell_porosity",
        # "fracture_porosity",
        # "cell_permeability",
        # "fracture_permeability",
        # "cell_thermal_conductivity",
        # "fracture_thermal_conductivity",
        "Psat",
        "Tsat",
        "liquid_molar_density",
        "liquid_molar_enthalpy",
        "liquid_dynamic_viscosity",
        "get_fluid_properties",
        "molar_density",
        "molar_enthalpy",
        "dynamic_viscosity",
        "is_context_locked",
        "lock_context",
        "unlock_context",
        "build_state",
        "all_states",
        # wells
        "nb_producers",
        "nb_injectors",
        "producers_information",
        "injectors_information",
        "producers_data",
        "injectors_data",
        "producer_perforations",
        "injector_perforations",
        "producers_perforations",
        "injectors_perforations",
        "NumNodebyProc",
        "NumFracbyProc",
        "NumWellInjbyProc",
        "NumWellProdbyProc",
        "retrieve_jacobian",
        "retrieve_right_hand_side",
        "retrieve_partitioning",
        "set_AMPI_cpp",
        "get_AMPI_nnz_cpp",
        "set_RHS_cpp",
    ]
)

_a_verifier = [
    # autonome
    "dump_array_in_fortran",
    "increment_first_column_in_fortran",
    # à trier
    "build_grid",
    "check_IncCV",
    "clear_all_neumann_contributions",
    "neumann_conditions",
    "dump_incv_info",
    "model_number_of_phases",
    "model_number_of_components",
    "model_number_of_contexts",
    "init_warmup",
    "init_phase2_partition",
    "init_phase2_build_local_mesh",
    "init_phase2_setup_contexts",
    "init_phase2_setup_VAG",
    "init_phase2_setup_solvers",
    "Residu_update_accumulation",
    "SolvePetsc_SetUp",
    "SolvePetsc_Init",
    "SolvePetsc_KspSolveIterationNumber",
    "SolvePetsc_dump_system",
    "SolvePetsc_Ksp_configuration",
    "SolvePetsc_Ksp_iterations",
    "SolvePetsc_ksp_solve",
    "SolvePetsc_check_solution",
    "SyncPetsc_GetSolNodeFracWell",
    "SyncPetsc_global_matrix_size",
    "SyncPetsc_local_matrix_size",
    "part_info",
    "SyncPetsc_rowcolnum",
    "SyncPetsc_colnum",
    "IncCV_SaveIncPreviousTimeStep",
    "IncCV_LoadIncPreviousTimeStep",
    "IncCVWells_PressureDrop",
    "DirichletContribution_update",
    "IncPrimSecd_update_secondary_dependencies",
    "LoisThermoHydro_compute",
    "Flux_DarcyFlux_Cell",
    "Flux_DarcyFlux_Frac",
    "Flux_FourierFlux_Cell",
    "Flux_FourierFlux_Frac",
    "Residu_compute",
    "Residu_reset_history",
    "Residu_RelativeNorm_local_closure",
    "Jacobian_ComputeJacSm",
    "NN_flash_all_control_volumes",
    "DefFlashWells_TimeFlash",
    "IncCVReservoir_NewtonRelax",
    "IncCV_NewtonIncrement",
    "IncPrimSecd_PrimToSecd",
    "Jacobian_GetSolCell",
    "gas_molar_enthalpy",
    "gas_molar_density",
    "gas_dynamic_viscosity",
]
