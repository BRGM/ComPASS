.. meta::
    :scope: version5

Setting up physical properties
==============================

As of now some properties can be defined on a mesh element or control volume basis,
whereas other physical properties are defined globally for the whole mesh.

.. The bindings are missing, no way to get the properties !
.. You can retrieve global properties using the corresponding `get_*` functions, e.g.:

.. .. code-block:: python

..     print("Rock volumetric heat capacity:", simulation.get_rock_volumetric_heat_capacity())


Gravity
-------

By default the gravity is 9.81. To modify it,
use :code:`model.gravity = the_gravity`.


Regionalized data
-----------------

Petrophysics
............

Properties such as permeability, porosity, thermal conductivity are defined
through the :code:`model.data` object.
It is possible to define the data using a constant value everywhere,
or setting a value in a zone, or giving an array...

.. code-block:: python

    model = Coats("physics_name")
    geom = GeomProperties(regular_mesh((Nx, Ny, Nz), extent=(Lx, Ly, Lz)))
    data = model.new_data(geom.mesh.dimension)

    data.petrophysics.cell_porosity[...] = 0.15
    data.petrophysics.cell_thermal_conductivity[...] = 3.0  # W/m/K
    # init the permeability everywhere
    data.petrophysics.cell_permeability[...] = 1e-12  # m^2
    # modify permeability values depending on the cell depth
    z_cell = Data(geom.cell_centers)[:, -1]
    data.petrophysics.cell_permeability[z_cell > H] = 1e-18


Rocktypes
.........

The rocktype is necessary to identify which
:ref:`capillary pressure and relative permeabilities
<Capillary pressure and relative permeabilities>`
are adapted at a certain location of the mesh.

All the cells must be initialized, the sites (cells and nodes in VAG)
can be deduced from cells value. Sites value can also be
specified if necessary.

.. code-block:: python

    # init with 0 everywhere, then modify some of them
    data.rocktypes[...] = 0
    z_cell = Data(geom.cell_centers)[:, -1]
    # using a condition on the depth of the cells centers
    data.rocktypes[z_cell > H] = 1

.. # or giving an array
.. all_cells = Zone(geom.mesh.cells)
.. assert len(all_cells) == len(rdomain)
.. data.rocktypes[all_cells] = rdomain



Fluid properties
----------------

Each physics comes with default physical properties,
they are listed :ref:`here <Available physics>`.

All the fluid properties can be adapted:
  - the molar density and the component molar mass, the volumetric
    mass density is updated with
    :math:`\rho^\alpha = \zeta^\alpha \sum\limits_{i\in\mathcal{C}} C_i^\alpha  M_i`,
  - the viscosity,
  - the molar enthalpy,
  - the saturation pressure.

You can choose one of the property already defined (be carreful that the
phase.s and component.s are compatible) or you can define your own function.
The property must be defined **with and without** the partial derivatives
**for each phase**.

When setting the property some checks are done by default to test that the
functions with and without the partial derivatives return the same value and
that the partial derivatives are correct. To desactivate the verifications, use the
:code:`check_derivatives=False` option.

Molar density, viscosity or molar enthalpy
..........................................

The property function without the partial derivatives returns the property
value; its only argument is the **phase state** :code:`X`.
The property function with the partial derivatives returns the property
value and has two arguments :
  - :code:`X` which contains the **phase state** *X.temperature*,
    *X.pressure* and *X.molar_fractions*,
  - :code:`dfdX` which is filled with the partial derivatives *dfdX.temperature*,
    *dfdX.pressure* and *dfdX.molar_fractions*.


.. code-block:: python

    def first_phase_property_with_derivatives(X, dfdX):
        dfdX.pressure = 1.0
        dfdX.molar_fractions.fill(0)
        dfdX.temperature = 0.3
        return 0.3 * X.temperature + X.pressure

    def first_phase_property_without_derivatives(X):
        return 0.3 * X.temperature + X.pressure

    first_phase_property = PhaseProperty(
        with_derivatives=first_phase_property_with_derivatives,
        without_derivatives=first_phase_property_without_derivatives,
    )


Then the phase properties must be set (follows an example with a two-phase physics)
where *property* is **molar_density**, **viscosity** or **enthalpy**:

.. code-block:: python

    model = Coats("diphasic")
    set_{property}_functions(
        model.properties.fluid,
        [
            first_phase_property,
            second_phase_property,
        ],
    )

The modification of the molar density will update by default the volumetric
mass density. Use the :code:`_update_volumetric_mass_density_functions`
optional argument to deactivate the automatic update.
It is useful for example if the components molar masses are also
modified, which also updates the densities.

.. code-block:: python

    set_molar_density_functions(
        model.properties.fluid,
        [
            first_phase_property,
            second_phase_property,
        ],
        _update_volumetric_mass_density_functions=False
    )

Component molar mass
....................

Use :code:`set_components_molar_mass` to modify the components molar masses,
follows an example with two components:

.. code-block:: python

    molar_masses = [29.0e-3, 0.018016]
    set_components_molar_mass(model.properties.fluid, molar_masses)

The volumetric mass density (expressed in kg/m^3) is set as
:math:`\rho^\alpha = \zeta^\alpha \sum\limits_{i\in\mathcal{C}} C_i^\alpha  M_i`.

.. warning::

    For implementation reasons, up to now you must modify the molar densities
    **before** updating the component molar masses.


Saturation pressure
...................

The functions of the saturation pressure take
a temperature state :code:`X` with only
the *X.temperature* value.

.. code-block:: python

    def psat_with_derivatives(X, dfdX):
        dfdX.temperature = 2 * X.temperature
        return X.temperature**2.0

    def psat_without_derivatives(X):
        return X.temperature**2.0

    psat = PhaseProperty(
        with_derivatives=psat_with_derivatives,
        without_derivatives=psat_without_derivatives,
    )

    set_psat_function(
        model.properties.fluid,
        psat,
        # to desactivate coherencies and derivatives validations
        check_derivatives=False,
    )

Rock properties
---------------

Rock volumetric heat capacity
.............................

For the time-being the rock volumetric heat capacity is defined
on a global basis (constant for the whole mesh).

The modification is done as follows:

.. code-block:: python

    model = Coats("diphasic")
    rho_rock = 2000 # kg/m^3 rock specific mass
    cp_rock = 800  # J/K/kg specific heat capacity
    set_rock_volumetric_heat_capacity(model.properties.rock, rho_rock * rho_cp)


Capillary pressure and relative permeabilities
...............................................

Particular regionalized properties are the capillary pressure
and the relative permeabilities
when there is at least two phases.
It is possible de define your own properties or to load
some already implemented.
It is regionalized via the :ref:`rocktype <Rocktypes>` defined for each
cell.

By default the capillary pressure is null, the phase relative permeability
is the square of the phase saturation :math:`kr^\alpha = (S^\alpha)^2`.

Use
:code:`physicalprop.set_relative_permeability.set_relative_permeability_functions`
and :code:`physicalprop.set_capillary_pressure.set_capillary_pressure_functions`
to set other capillary pressure or relative permeabilities.

Some utilities are available,
for example you can use the already implemented van Genuchten
capillary pressure and relative permeabilities functions:

.. code-block:: python

    from physicalprop.set_van_Genuchten_rock_properties import (
        add_van_Genuchten_rock_properties,
    )
    model = Coats("the_physics")
    # define capillary pressure and relative permeability
    # with Van Genuchten formula for the domain where rocktype=1
    add_van_Genuchten_rock_properties(
        model.properties.rock, rocktype=1, Pr=15.0e6, Slr=0.4, Sgr=0, n=1.49
    )


Developpement utilities
-----------------------

Test python properties
......................

To test the implementation of the property function, you can use the
:code:`PhaseStateStruct` class whose constructor takes the number of components.

.. code-block:: python

    from physics.physical_state import PhaseStateStruct
    p = 1.0 * bar
    T = 280.0  # K
    C = np.array([0.0, 1.0])
    phase_state_type = PhaseStateStruct(nb_components)
    # creates X such that X.pressure=p, X.temperature=T and X.molar_fractions=C
    X = phase_state_type.Xalpha(p, T, C)
    # creates the good shape object for partial derivatives
    dfdX = phase_state_type.empty_Xalpha()
    property_value = first_phase_property_with_derivatives(X, dfdX)
    print("The property is equal to ", property_value)
    print("The derivative with respect to the pressure is ", dfdX.pressure)


.. The file :download:`call_python_viscosity.py <../test/unit/call_python_viscosity.py>`
.. presents more examples about how to call the property function.

Constant fluid properties
.........................

If a fluid property is constant, use :code:`constant_fluid_physical_property`.
It initialized the functions with and without the derivatives.

.. code-block:: python

    from physicalprop.utils import constant_fluid_physical_property

    simulation.set_viscosity_functions(
        property_functions=[
            constant_fluid_physical_property(2.0e-5),
            constant_fluid_physical_property(1.0e-3),
        ]
    )
