.. meta::
    :scope: version5

Import or create a mesh
=======================

The present section gives information on how to generate or import a mesh in ComPASS.
It is currently possible to generate 3D cartesian meshes using the
`icmesh <https://gitlab.com/compass/compass-v5/icmesh>`_ software brick
which creates the mesh and the
`geom-traits <https://gitlab.com/compass/compass-v5/geom-traits>`_ brick
which completes the mesh with geometry utilities.


Creating a regular Cartesian Mesh
---------------------------------

Create a cartesian mesh with uniform cells size
in the *x*, *y*, and *z* directions.

.. code-block:: python

    from icmesh import regular_mesh
    Nx = 1; Ny = 5; Nz = 64  # number of cells
    Lx = 1.0; Ly = 3.0; Lz = 320.0  # cell lengths (in m)
    Ox = 0.0; Oy = -3.0; Oz = 10.0  # origin position (in m)

    # origin is optional: (0,0,0) by default
    the_mesh = regular_mesh((Nx, Ny, Nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz))


Creating a Scottish Cartesian Mesh
----------------------------------

Create a cartesian mesh with non-constant cells size.
:code:`icmesh.scottish_mesh` takes the cell lengths array
in each direction.

.. code-block:: python

    from icmesh import scottish_mesh

    def scottish_step(N, L):
        factor = 2 * L / N / (N + 1)
        return np.array([factor * (i + 1) for i in range(N)])

    Nx = 1; Ny = 5; Nz = 64  # number of cells
    Lx = 1.0; Ly = 3.0; Lz = 320.0  # cell lengths (in m)
    Ox = 0.0; Oy = -3.0; Oz = 10.0  # origin position (in m)

    # origin is optional: (0,0,0) by default
    the_mesh = scottish_mesh(
        scottish_step(Nx, Lx), scottish_step(Ny, Ly), scottish_step(Nz, Lz),
        origin=(Ox, Oy, Oz),
    )


Visualize the mesh
------------------

The visualization of the mesh in Paraview is facilitated using
:code:`icmesh.dump.dump_vtu_3d(the_mesh, "visu_file.vtu")`
which creates the vtu corresponding file.
You can add values carried by the mesh elements using *pointdata*,
*celldata* and/or *fielddata* such as:

.. code-block:: python

    from icmesh.dump import dump_vtu_3d
    assert len(rdomain) == len(the_mesh.cells)
    dump_vtu_3d(the_mesh, "visu_file.vtu", celldata={"domain": rdomain})


Add geometry utilities
----------------------

When loading a mesh for a whole simulation, it is necessary to add
geometry considerations using the
`geom-traits <https://gitlab.com/compass/compass-v5/geom-traits>`_ compass
software:

.. code-block:: python

    from icmesh import regular_mesh
    from geom_traits import GeomProperties
    Nx = 1; Ny = 5; Nz = 64  # number of cells
    Lx = 1.0; Ly = 3.0; Lz = 320.0  # cell lengths (in m)

    geom = GeomProperties(regular_mesh((Nx, Ny, Nz), extent=(Lx, Ly, Lz)))
    # the mesh is accessible through geom.mesh
    dump_vtu_3d(geom.mesh, "visu_file.vtu")
    # access mesh geometry utilities, for example
    cell_centers = geom.cell_centers
    cell_measures = geom.cell_measures
