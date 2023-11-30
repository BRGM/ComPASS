Tutorial
========

The following is a tutorial starting from a very simple script
on which we will progressively add more complex physical processes.

Step 1: Execution and visualization
-----------------------------------

Let's start with the following script:

.. literalinclude:: ../tutorial/water2ph_vertical_column.py
   :language: python
   :linenos:

Download file:
:download:`water2ph_vertical_column.py <../tutorial/water2ph_vertical_column.py>`


Execution:

.. code-block:: shell

  cd ./tutorial/
  python3 water2ph_vertical_column.py


* What contains this simulation (mesh, petrophysics, boundary condition, initial states, ...)?

* What are the output in the terminal?

* To visualize using Paraview:

 Add :code:`simulation.postprocess()` at the end of the script and execute again.

 Or execute in the terminal:

.. code-block:: shell

  python3 -m ComPASS.postprocess -s -C output-water2ph_vertical_column/

The :func:`postprocess` function creates the states.pvd file in the
output-water2ph_vertical_column/paraview/ directory. You can precise
an other directory name in :code:`ComPASS.set_output_directory_and_logfile()`.

* By default, information about the simulation and the convergence
  are stored for each time step in the
  output-water2ph_vertical_column/time_step_log directory
  and a summary in output-water2ph_vertical_column/timeloop_log.yaml.
  Some scripts to analyse the information are in the ../test/utilities/
  directory.

For exemple you can execute:

.. code-block:: shell

  python3 ../test/utilities/timeloops_analysis.py ./output-water2ph_vertical_column

Step 2: Improve the script
--------------------------

* Remove the :code:`set_initial_states` function to use the build-in
  method :code:`simulation.build_state()`
  An example is in the :ref:`setting-initial-values` section. Don't forget to set the Dirichlet nodes,
  they are not contained in :code:`simulation.all_states()`. Indeed, :code:`simulation.all_states()`
  contains all the states needed to initialize the VAG sites : the nodes, the cells
  and the fracture faces if any, but it does not depends on the boundary conditions.

* Transform the bottom boundary into a :ref:`Neumann BC<neumann_faces_bc>` with
  heat flux only (bottom_heat_flux = :math:`0.08  W.m^{-2}`).


Step 3: Change the Physics
----------------------------

* Change the physics from *water2ph* to *diphasic*. A description of
  the disponible physics is in the :ref:`documentation<physics_section>` section.

How many phases and components do we use now? Which context exist?
Which differences in the output in the terminal compared to previously?

* Add gravity changing the value in :code:`simulation.set_gravity(g)`.

* Initialize the domain with hydrostatic pressure :
  p = p0 + rho * gravity * (zTop - z), with rho = 1e3.
  A very similar example is done in the :ref:`setting-initial-values` section.

* Change the top boundary condition to impose a Dirichlet diphasic state with
  Sg = 0.5 at p0 and T0. You can add :code:`Sg=0.5` in :code:`simulation.build_state`.
  Is the air component present in the simulation? Why? You can print the state as follows:

.. code-block:: python

  Xdir = build_state(...)
  print(Xdir)

Step 4: Van Genuchten capillary pressure and relative permeabilities
--------------------------------------------------------------------

* Plot the van Genuchten Pc and the kr by executing the capillary file
  and the relative permeabilities file.
  The already implemented capillary pressures are in the directory
  ./ComPASS/petrophysics/models,
  the relative permeabilities are in ./tutorial/data/van_genuchten_kr.py

* Add a van Genuchten capillary pressure and the corresponding relative permeabilities
  (some :ref:`documentation<pc_kr>`).
  To avoid lots of time step failures,
  you might change the initial_timestep and the increase_factor.

  What is the impact on the gas saturation?

.. warning::

  Pc and kr should be set after the partition (done in :code:`simulation.init`)
  and before the :func:`build_state` function of the diphasic physics
  (because it calls the phase pressure function).


* Set the cell :ref:`rocktypes<setting_rocktypes>` : the top half of the mesh has a rocktype = 1 (COX),
  elsewhere rocktype = 2 (CCT) (these two rocktypes values are
  implemented in the van Genuchten Pc and kr).
  You may need the :func:`simulation.compute_global_cell_centers()` method to get the
  coordinates of all the cells centers.

  What is the impact on the gas saturation?


Step 5: Add a vertical fracture
-------------------------------

* I strongly advise (in a first step) to deactivate the capillary pressure and
  the relative permeabilities to help the convergence.

* Change the mesh to have Lx = 400, Ly = 100, H = 800.
  Also change the discretization with nx = 20, ny = 4, nz = 40.
  The number of cells has increased, thus you might prefer to run in parallel:

.. code-block:: shell

    mpirun -n 2 python3 water2ph_vertical_column.py

* Add a vertical fracture at (y==0) & (z < - H/2) of thickness 0.1 m,
  with permeability = 1e-12 m^2, porosity = 0.5
  and thermal conductivity = 2 W/m/K. Have a look at the :ref:`fracture
  documentation<fractures_sec>`, you may need
  :code:`simulation.compute_global_face_centers()` to obtain the coordinates
  of all the face centers.

* Change the Neumann boundary condition at the bottom to apply a
  :ref:`Neumann flux at the fracture edges<frac_edges_bc>` only.
  Apply a molar flux of :math:`0.1 mol.m^{-2}.s^{-1}` on the water component
  (the water component is the second one) and the corresponding heat flux
  using the :func:`simulation.liquid_molar_enthalpy()` function. In the diphasic
  physics this function takes the following arguments: the liquid pressure, the temperature and the phases
  molar fraction (use :code:`pure_phase_molar_fraction = [[0, 1], [1, 0]]`).
  The heat flux is equal to the molar flux times the liquid molar enthalpy.

  Does the Newton algorithm converge? Why?

* To obtain the convegence, you need to change the maximum number of iterations of the
  :ref:`Newton algorithm<setting_newton>`
  using the :func:`Newton` function
  from :code:`ComPASS.newton`.

* You can activate again the capillary pressure and
  the relative permeabilities.


Step 6: Atmospheric boundary condition
--------------------------------------

* Deactivate the capillary pressure and
  the relative permeabilities to help the convergence.

* Change the top boundary condition: the top nodes with coordinates x <= 0
  remains Dirichlet BC (with Sg=0.5, p=p0, T=T0),
  the top faces with face center coordinates x >= 0 becomes :ref:`atmospheric BC<atmBC>`.
  Initialize the porous media **nodes** where the atmospheric BC is imposed
  with Sg=0.5, p=p0 and T=T0,
  it corresponds to the :code:`simulation.Context.diphasic_FF_no_liq_outflow`
  context.

* Modify the far-field value of the
  :ref:`atmospheric boundary condition<far_field_atmBC>`
  to set the temperature to 20C.

* What happends if you add the capillary pressure and
  the relative permeabilities?

* The convergence is difficult to obtain, in particular the first time step.
  Then a solution is to call the
  :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>` function
  a first time without the van Genuchten Pc nor kr over a very small
  number of time steps, to do so add :code:`nitermax = 5` in the argument
  list of :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>`.
  It will stop the time loop after 5 iterations, which will compute
  an initialization.
  After the call to the
  :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>` function,
  set the Pc and kr, don't forget to initialize again the Dirichlet nodes
  (the :func:`simulation.build_state` function depends on the Pc)
  and call again the
  :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>` function
  to simulate the total time.

  Remark: the output of the :func:`simulation.standard_loop` is the time
  reached. It is useful sometimes so set that time as the initial time
  of the next call of the time loop.

.. code-block:: python

  current_time = simulation.standard_loop(...)
  print(current_time)
  simulation.standard_loop(initial_time=current_time, ...)


* In complex execution it is often necessary to initialize the simulation
  with a specific state to reach an equilibrium, then to adapt or change
  some parameters and run the complex execution.
  In such a case, you can also :ref:`reload any previous simulation state
  <setting-initial-values>`
  from an output directory using the :func:`simulation.reload_snapshot` method.
