.. meta::
    :scope: version5

Main structure of the scripts
=============================

A simulation script is composed of several main parts.

Script is in `test_diphasic_Dirichlet.py. <https://gitlab.com/compass/compass-v5/compass-coats/-/blob/main/test/test_diphasic_Dirichlet.py?ref_type=heads>`_

* Creating or loading a :ref:`mesh <Import or create a mesh>` together with
  geometry utilities.

  .. code-block:: python

      from icmesh import regular_mesh
      from geom_traits import GeomProperties
      Nx = 1; Ny = 5; Nz = 64  # number of cells
      Lx = 1.0; Ly = 3.0; Lz = 320.0  # cell lengths (in m)
      geom = GeomProperties(regular_mesh((Nx, Ny, Nz), extent=(Lx, Ly, Lz)))
.. A mesh can be created explicitely or loaded from a file.

* Choosing the :ref:`numerical scheme <Numerical scheme>`.

  .. code-block:: python

    from compass_coats.schemes import TPFA, VAG
    scheme_def = TPFA()  # or VAG()

* Loading a model from one of the :ref:`available physics <Load the physics>`, e.g.:

  .. code-block:: python

    from compass_coats.models import Coats
    model = Coats("diphasic")

  To change the gravity :code:`model.gravity = the_gravity`.

* Initializing the :ref:`simulation properties <Regionalized data>`
  (petrophysics, rocktypes...):

  .. code-block:: python

    data = model.new_data(geom.mesh.dimension)
    data.petrophysics.cell_porosity[...] = 0.15
    data.petrophysics.cell_permeability[...] = 1e-12
    data.petrophysics.cell_thermal_conductivity[...] = 3.0
    # rocktypes need to be initialized over cells
    data.rocktypes[...] = 0

* Setting-up :ref:`initial values<Setting up initial values>`
  and :ref:`boundary conditions<Setting up boundary conditions>`.

  .. code-block:: python

    # initial state and BC
    init_ctx, Xinit = model.utils.equilibriate_state(
        context="liquid",
        p=10 * bar,
        T=degC2K(15),
        Cal=1e-4,
    )
    top_ctx, Xtop = model.utils.equilibriate_state(
        context="diphasic",
        p=11 * bar,
        T=degC2K(35),
        Sg=0.1,
    )
    # init all states with initial values
    data.initial_states[...] = Xinit
    data.initial_contexts[...] = init_ctx
    top = top_boundary(geom)
    data.boundary_conditions.Dirichlet_states[top] = Xtop
    data.boundary_conditions.Dirichlet_contexts[top] = top_ctx

* Solving the :ref:`temporal problem <Time loop execution>` using the
  :code:`compass-coats.Standard_time_loop` class.

  .. code-block:: python

    from compass_coats.output_visu import output_directory
    from compass_coats.standard_time_loop import Standard_time_loop

    visu_dir = output_directory(__file__)
    time_loop = Standard_time_loop(
        geom=geom,
        model=model,
        scheme=scheme_def,
        data=data,
        output_dir=visu_dir,
    )
    # if necessary, adapt the Newton or timestepper coefficients
    time_loop.loop.timestepper.step_solver.maxiter = 25
    time_loop.loop.timestep_manager.increase_factor = 1.5

    solution, tick = time_loop.run(
        initial_step=1e4,
        final_time=100 * day,
        output_every=10,
    )

* Simulation results can be prostprocess with compass_coats.postprocess module
  and the underlying :func:`~compass_coats.postprocess.postprocess` function.

  .. code-block:: python

    from compass_coats.postprocess import postprocess
    postprocess(visu_dir, time_unit="day")
