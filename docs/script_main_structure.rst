Main structure of the scripts
=============================

A simulation script is composed of 4 main parts:

1. Defining a simulation from one of the available physics, e.g.:
:code:`simulation=ComPASS.load_eos("water2ph")`

and initializaing variables, including creating and loading a global mesh.
A mesh can be created explicitely or loaded from a file.
Several io functions are available from the MeshTools package and a MeshTools object can be
passed to the :code:`simulation.init` method (cf. phase 2).

2. Definition of fractures, simulation properties (permeabilities), partitioning of the mesh.
All of this is done calling the :code:`simulation.init` method
(cf. :func:`~ComPASS.simulation.base.init` documentation).
After step 2 the mesh is partitioned and the global mesh is no longer available.
Each proc works on its own subdomain.

3. Setting-up initial values and physical states (cf. various examples in 
:ref:`Setting-up initial conditions`_)

4. Solving the temporal problem, i.e. making or using one of the provided time loops.
Cf. for example the :func:`ComPASS.timeloops.standard_loop` function.


Simulation results can be prostprocess with ComPASS.prostprocess module.
If you use docker container you can also use the :code:`posprocess` subcommand.
