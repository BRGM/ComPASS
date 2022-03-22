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


* What contains this simulation (mesh, petrophysics, boundary condition, initial states, ...) ?

* What are the output in the terminal ?

* Visualization:

 Add :code:`simulation.postprocess()` at the end of the script and execute again.

 Or execute in the terminal:

.. code-block:: shell

  python3 -m ComPASS.postprocess -s -C output-water2ph_vertical_column/

The :func:`postprocess` function creates the states.pvd file in the
output-water2ph_vertical_column/paraview/ directory.



Step 2: Improve the script
--------------------------

* Remove the :code:`set_initial_states` function to use the build-in
  method :code:`simulation.build_state()`
  An example is in the :ref:`setting-initial-values` section. Don't forget to set the Dirichlet nodes,
  they are not contained in :code:`simulation.all_states()`.

* Transform the bottom boundary into a Neumann BC with
  heat flux only (bottom_heat_flux = 0.08 J/m^2/s).


Step 3: Change the EOS
----------------------

* Change the Equations of states from water2ph to diphasic.

How many phases and components do we use now ?

* Add gravity.

* Initialize the domain with hydrostatic pressure :
  p = p0 + rho * gravity * (zTop - z), with rho = 1e3

* Change the top boundary condition to impose a diphasic state (Sg = 0.5) at p0 and T0.
 Is the air component present in the simulation ? Why ?


Step 4: Van Genuchten capillary pressure and relative permeabilities
--------------------------------------------------------------------

* Plot the Van Genuchten Pc and the kr (already coded in the source code).

* Add a van Genuchten capillary pressure and the corresponding relative permeabilities.
  To avoid lots of time step failures,
  you might change the initial_timestep and the increase_factor.

  What is the impact on the gas saturation ?

.. warning::

  Pc and kr should be set after the partition (done in :code:`simulation.init`)
  and before the :func:`build_state` function of the diphasic eos
  (because it calls the phase pressure function).


* Set the cell rocktypes : the top half of the mesh has a rocktype = 1 (COX),
  elsewhere rocktype = 2 (CCT) (these two rocktypes values are implemented in the van Genuchten Pc).

  What is the impact on the gas saturation ?


Step 5: Change the mesh and add a vertical fracture
---------------------------------------------------

* I advise (in a first step) to deactivate the capillary pressure and
  the relative permeabilities to help the convergence.

* Change the mesh to have Lx = 400, Ly = 100, H = 800.
  Also change the discretization with nx = 20, ny = 4, nz = 40.
  The number of cells has increased, thus you might prefer to run in parallel:

.. code-block:: shell

    mpirun -n 2 python3 water2ph_vertical_column.py

* Add a vertical fracture at (y==0) & (z < - H/2) of thickness 0.1 m,
  with permeability = 1e-12 m^2, porosity = 0.5
  and thermal conductivity = 2 W/m/K.

* Change the Neumann boundary condition at the bottom to apply a flux
  at the fracture edges only.
  Apply a molar flux of 0.1 mol/m^2/s on the water component
  (the water component is the second one) and the corresponding heat flux
  using the :func:`simulation.liquid_molar_enthalpy` function.

  Does it converge ? Why ?

* Change the maximum number of Newton iterations using the :func:`Newton` function
  from :code:`ComPASS.newton`;
  and the :code:`newton` keyword in :func:`simulation.init`.

* You can activate again the capillary pressure and
  the relative permeabilities.


Step 6: Atmospheric boundary condition
---------------------------------------------------

* Deactivate the capillary pressure and
  the relative permeabilities to help the convergence.

* Change the top boundary condition: the top nodes with coordinates x <= 0
  remains Dirichlet BC (with Sg = 0.5, p=p0, T=T0),
  the top faces with face center coordinates x >= 0 becomes atmospheric BC.
  Initialize the porous media **nodes** where the atmospheric BC is imposed
  with Sg=0.5, p=p0 and T=T0,
  it corresponds to the :code:`simulation.Context.diphasic_FF_no_liq_outflow`
  context.

.. warning::

  Up to now the Dirichlet node states must be set after the node states
  where the atmospheric BC is imposed.

* Modify the far-field value of the atmospheric boundary condition
  to set the temperature to 20C.

* What happends if you add the capillary pressure and
  the relative permeabilities ?

* To help the convergence, call the :func:`simulation.standard_loop`
  function before adding the Pc and kr, with a small maximum number of
  iterartions (for example :code:`nitermax = 5`), it will initialize
  the simulation.
  Set the Pc and kr after the first :func:`simulation.standard_loop` call,
  initialize again the Dirichlet nodes
  (the :func:`simulation.build_state` function depends on the Pc)
  and call again the :func:`simulation.standard_loop` function to simulate
  the total time.
