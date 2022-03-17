Setting up physical properties
==============================

As of now some properties can be defined on a mesh element or control volume basis,
whereas other physical properties are defined globally for the whole mesh.

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
    )


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
