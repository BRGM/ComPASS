.. _setting_physical_properties:

Setting up physical properties
==============================

As of now some properties can be defined on a mesh element or control volume basis,
whereas other physical properties are defined globally for the whole mesh.


Global properties
-----------------

For the time-being several properties are defined on a global basis
calling global functions for the whole mesh.

Global properties are the following:
  - gravity
  - fracture thickness
  - rock volumetric heat capacity


For example setting up the rock volumetric heat capacity can be done as follows:


.. code-block:: python

    rho_rock = 2000 # kg/m^3 rock specific mass
    cp_rock = 800  # J/K/kg specific heat capacity
    simulation.set_rock_volumetric_heat_capacity(rho_rock * rho_cp)

You can retrieve global properties using the corresponding `get_*` functions, e.g.:

.. code-block:: python

    print("Rock volumetric heat capacity:", simulation.get_rock_volumetric_heat_capacity())

.. warning::

    The global functions to set global properties are to be called before
    the :func:`simulation.init <ComPASS.simulation.base.init>` method
    that will distribute all properties to the partitioned mesh.

Regionalized properties
-----------------------

Properties such as permeability, porosity, thermal conductivity are defined
through the
:func:`simulation.init <ComPASS.simulation.base.init>` function and the corresponding keywords for
cell or fracture elements.


For example setting up the reservoir thermal conductivity can be done as follows:

.. code-block:: python

    K_reservoir = 2 # bulk thermal conductivity in W/m/K
    simulation.init(
        ...,
        cell_thermal_conductivity=K_reservoir,
        ...,
    )

Mesh objects identification after the distribution
-------------------------------------------------

With the parallelism, the indexes of the mesh objects
(cells, nodes, faces) change after the distribution of the mesh
done in the
:func:`simulation.init <ComPASS.simulation.base.init>` function.
ComPASS includes a tool named *flags* to identify the objects
after the distribution. It is usefull for example when you tag
objects using a mesh generator and you want to use the information
after the distribution. The flags are not used elsewhere in the ComPASS
code, it is for you to track some mesh information.

The first step is to set the global flags using the
:code:`set_global_flags` keyword in the
:func:`simulation.init <ComPASS.simulation.base.init>` function. The flags
contain one integer by object, which is initialized to 0.
Then you can retrieve the flags with the local indexes of the mesh objects.
For example:

.. code-block:: python

    gallery_flag = 3
    def set_flags():
        nodeflags = simulation.global_nodeflags() #  already init with 0
        vertices = simulation.global_vertices() #  all vertices coordinates
        gallery_vertices = on_zmin(grid)(vertices) #  bool array
        nodeflags[gallery_vertices] = gallery_flag

    simulation.init(
        ...,
        set_global_flags=set_flags,
        ...,
    )

    nodeflags = simulation.nodeflags()
    Xgal = simulation.build_state(simulation.Context.gas, p=pgal, T=Tgal)
    simulation.node_states().set(nodeflags == gallery_flag, Xgal)


Capillary pressure and relative permeabilities
----------------------------------------------

Particular regionalized properties are the capillary pressure
and the relative permeabilities when there is at least
two phases.
It is possible de define your own laws or to load
one already implemented.
It is regionalized via the rocktype defined for each
cell and fracture face. Then the rocktype is a variable
used by the code that you can initialize (by default is 1).

.. _setting_rocktypes:

Rocktypes
.........

The rocktype is set using the
:code:`set_global_rocktype` keyword in the
:func:`simulation.init <ComPASS.simulation.base.init>` function.
For example:

.. code-block:: python

    def select_global_rocktype():
        # you can define the rocktype (for example depending on the geometry)
        cell_centers = simulation.compute_global_cell_centers()
        COX = cell_centers[:, 1] > Lx / 2  # right half
        CCT = cell_centers[:, 1] <= Lx / 2.0
        cellrocktype = simulation.global_cell_rocktypes()
        cellrocktype[COX] = 1
        cellrocktype[CCT] = 2

        # or you can rely on values you already initialized in the flags
        # Careful : the size of global_fracture_rocktypes is NbFace !
        faceflags = simulation.global_faceflags()
        fracrocktype = simulation.global_fracture_rocktypes()
        fracrocktype[:] = np.stack((faceflags, faceflags), axis=-1)

    simulation.init(
        ...,
        set_global_rocktype=select_global_rocktype,
        ...,
    )

For non-isothermal simulation,
:code:`global_cell_rocktypes` (also :code:`global_fracture_rocktypes`) is
composed of two values, the first one for the capillary pressure, the seconde
one for a thermal use.
By default the rocktype values are 1.

.. warning::
    All the faces are in :code:`global_fracture_rocktypes`, not only the
    fracture faces.

.. _pc_kr:

Capillary pressure
..................

By default the capillary pressure is null.
The capillary functions are in the ComPASS/petrophysics/models directory and you can
define new ones in the python language.
For example, you can use the already implemented van Genuchten capillary
function as follows:

.. code-block:: python

    simulation.set_vanGenuchten_capillary_pressure()

Relative permeabilities
.......................

Using the rocktypes, it is also possible to define and regionalize
the relative permeability of each phase with the
:func:`simulation.set_kr_functions<ComPASS.simulation.base.set_kr_functions>`
function.

For example:

.. code-block:: python

    from data.van_genuchten_kr import kr_functions
    simulation.set_kr_functions(kr_functions)

By default (when the EOS contains two phases)
the relative permeability of each phase is the
square of the phase saturation.
