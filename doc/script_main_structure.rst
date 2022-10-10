Main structure of the scripts
=============================

A simulation script is composed of 4 main parts:

1. Defining a simulation from one of the available physics, e.g.:
:code:`simulation=ComPASS.load_physics("water2ph")`

and initializing variables, including creating and loading a global mesh.
A mesh can be created explicitely or loaded from a file.
Several io functions are available from the MeshTools package and a MeshTools object can be
passed to the :code:`simulation.init` method (cf. phase 2).

2. Definition of fractures, simulation properties (permeabilities), partitioning of the mesh.
All of this is done calling the :code:`simulation.init` method
(cf. :func:`~ComPASS.simulation.base.init` documentation).
After step 2 the mesh is partitioned and the global mesh is no longer available.
Each proc works on its own subdomain.

3. Setting-up initial values and physical states (cf. various examples in
:ref:`setting-initial-values`).

4. Solving the temporal problem, i.e. making or using one of the provided time loops.
Cf. for example the :func:`ComPASS.timeloops.standard_loop` function.


Simulation results can be prostprocess with ComPASS.prostprocess module
and the underlying :func:`~ComPASS.postprocess.postprocess` function.
If you use docker container you can also use the :code:`postprocess` subcommand.

Physical units
--------------
.. _units:

The
`International System of Units <https://en.wikipedia.org/wiki/International_System_of_Units>`_
is used throughout the code without specifying units.

For example:

   - distances are expressed in meters
   - permeabilities are expressed in :math:`m^2`
   - porosity has no units
   - times are expressed in seconds
   - ...

To help users, some precomputed quantities are available in
the `ComPASS.utils.units <https://github.com/BRGM/ComPASS/blob/v4.4.x/ComPASS/utils/units.py>`_ module.

For example, considering any function :code:`f(t)` that is expecting an argument in seconds,
just import another duration from :code:`ComPASS.utils.units` to call :code:`f`:

.. code-block :: python

  from ComPASS.utils.units import year

  f(10*year)

The examples in this section demonstrate how these quantities can be used in simulation scripts.
