.. meta::
    :scope: version4

Import or create a mesh
=======================

The present section gives information on how to generate or import a mesh in ComPASS.
The mesh is a key input of the :func:`simulation.init <ComPASS.simulation.base.init>` method.
It is possible to dump the global mesh before partitioning and properties distribution using
the keyword :code:`dump_mesh_before_distribution` which defaults to :code:`False`.


Creating a Cartesian Mesh
-------------------------

It is possible to create a cartesian uniform mesh using the :code:`Grid` from the :code:`ComPASS.Grid` module.
Many scripts uses this grid generator, for example :download:`linear_vertical_column.py <../test/bulk/linear_vertical_column.py>`.

Creating an extruded sector
---------------------------
The command :code:`extruded_sector` meshes an horizontal angular sector and extrudes it in 3D.
It is in the :code:`ComPASS.utils.angular_sector` module and it relies on MeshTools functions.

Example: :download:`MO1.py <../test/baseline/MO1/MO1.py>`


Using Salome
------------


Using MeshTools
---------------

Several io functions are available from the MeshTools package and a MeshTools object can be
passed as mesh to the :func:`simulation.init <ComPASS.simulation.base.init>` method.


Using Petrel/Eclipse grid format
--------------------------------

A faulted Petrel grid can be import and remeshed using CGAL under the hood using the :code:`PetrelGrid` from the :code:`ComPASS.utils.petrel` module.
The faces of the Petrel cells are split on their discontinuities to keep the number of cells constant (faster computation).
The permeability tensors are extracted from the Petrel file using the :code:`PetrelGrid.permeability` property.

More information and examples in the :ref:`example scripts section<eclipse_grid>`.

Fractures
---------
The fractures are objects of codimension 1, they are defined as a set of faces.
The geometry of the fractures is given with the keyword :code:`fracture_faces` in :func:`simulation.init <ComPASS.simulation.base.init>`.
It can be a list of face indexes, or a mask over the faces.

Wells
-----
Each well consists in a set of edges of the mesh. A list of well objects containing the geometry is passed via the keyword :code:`wells` in :func:`simulation.init <ComPASS.simulation.base.init>`.
More information in the :ref:`wells section<wells_introduction>`.
