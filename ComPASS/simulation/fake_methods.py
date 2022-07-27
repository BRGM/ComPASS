from .init import init

from ..timeloops import standard_loop

from ..timestep import make_one_timestep

from ..wells.wells import (
    create_well_from_segments,
    create_single_branch_well,
    create_vertical_well,
    get_well_data,
    get_well_perforations_state,
    get_wellhead,
    set_well_property,
    close_well,
    open_well,
    well_production_history,
    well_injection_history,
    close_perforations,
    set_well_model,
    check_well_geometry,
)

from ..utils.grid import (
    vertical_boundaries,
    bottom_boundary,
    top_boundary,
    all_boundaries,
)

from ..utils.phase_computations import total_phase_volume
from ..utils.various import (
    phases,
    components,
    contexts,
    states_locations,
    mass_fluxes_locations,
    enthalpy_fluxes_locations,
    reload_snapshot,
)
from ..utils.initialization import hydrostatic_pressure_profile

from .bc import (
    set_global_dirichlet_nodes,
    clear_dirichlet_nodes,
    set_dirichlet_nodes,
    reset_dirichlet_nodes_states,
    reset_dirichlet_nodes,
    set_freeflow_faces,
    reset_freeflow_faces,
    get_freeflow_nodes,
    get_freeflow_nodes_area,
)

from ..io.mesh import write_mesh, write_polyhedra_vtu_mesh

from .utils import facenodes, postprocess, eos_name, collect_all_edges

from ..petrophysics.kr import set_kr_functions
from ..petrophysics.phase_pressure import set_phase_pressure_functions
from ..petrophysics.capillary import set_liquid_capillary_pressure
from ..petrophysics.models.vanGenuchten import set_vanGenuchten_capillary_pressure

from ..linalg.factory import linear_solver
from ..newton import default_Newton
