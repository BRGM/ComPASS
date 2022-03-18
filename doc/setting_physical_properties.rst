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

Capillary pressure and relative permeability
--------------------------------------------

One particular regionalized property is the capillary pressure.
It is possible de define your own capillary pressure or to load
one already implemented.
The capillary pressure is regionalized via the rocktype defined for each
cell and fracture face. The rocktype is set using the
:code:`set_global_rocktype` keyword in the
:func:`simulation.init <ComPASS.simulation.base.init>` function.
For example:

.. code-block:: python

    def select_global_rocktype():
        cellflags = simulation.global_cellflags()
        cellrocktype = simulation.global_cell_rocktypes().reshape((-1, 2))
        cellrocktype[:] = np.stack((cellflags, cellflags), axis=-1)

        # Careful : the size of global_fracture_rocktypes is NbFace !
        faceflags = simulation.global_faceflags()
        fracrocktype = simulation.global_fracture_rocktypes().reshape((-1, 2))
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

And, for example, you can use the already implemented van Genuchten capillary
function as follows:

.. code-block:: python
    simulation.set_vanGenuchten_capillary_pressure()

The capillary functions are in the ComPASS/petrophyscs/models directory and you can
define new ones in the python language.
By default the capillary pressure is null.

Using the rocktypes, it is also possible to define and regionalize
the relative permeability of each phase with the
:func:`simulation.set_kr_functions<ComPASS.simulation.base.set_kr_functions>`
function.

For example:

.. code-block:: python

    from data.van_genuchten_kr import kr_functions
    simulation.set_kr_functions(kr_functions)
