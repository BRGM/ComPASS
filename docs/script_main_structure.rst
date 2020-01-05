
Main structure of the scripts
=============================

A simulation script is composed of 4 main parts:

1. Initialization of variables, including creating and loading a global mesh.

2. Definition of fractures, simulation properties (permeabilities), partitioning of the mesh.
All of this is done calling the :meth:`ComPASS.init` method.

3. Setting-up initial values and physical states.

4. Solving the temporal problem, i.e. making or using a time loop.

After step 2 the mesh is partitioned and the global mesh is no longer available.
Each proc works on its own subdomain.

Simulation results can be prostprocess with ComPASS.prostprocess module.


Phase 1
-------

A mesh can be created explicitely or loaded from a file.
Several io functions are available from the MeshTools package and a MeshTools object can be
passed to the :meth:`ComPASS.init` method (phase 2).

Phase 2
-------

All of its phase is done by calling the :meth:`ComPASS.init` method.
Most of the work will be performed on the master processor.
At the end of the method the mesh is partitioned and properties are distributed
to all procs.

Initialize many simulation properties and distribute the mesh.
Before a call to ComPASS.init, the mesh is global in the sense that
there is only one mesh on the master proc.
After the execution of ComPASS.init the mesh is distributed, i.e.
there is as many local meshes (with possibly ghost elements) as procs.


:param mesh: the mesh the simulation is run on, it can be a grid generated with ComPASS.Grid,
    or a MeshTools object. MeshTools provides several helpers modules to load and modify meshes.

:param wells: a python sequence of well objects
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
