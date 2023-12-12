.. meta::
    :scope: version5

Setting up initial values
=========================

Before any time-loop the numerical scheme sites values
must be initialized.

It is necessary to set both
:code:`data.initial_states` and :code:`data.initial_contexts`.
They are carried by mesh objects.

The initial values must be given over all the sites, ie all the **cells** must
be initialized with :ref:`TPFA<TPFA scheme>`;
all the **cells and nodes** must be initialized
with :ref:`VAG<VAG scheme>`.

.. Have a look at :ref:`this example<Setting up initial inputs>`
.. to discover different ways to initialize them.

.. code-block:: python

    data = model.new_data(geom.mesh.dimension)
    # build a liquid state at fixed pressure and temperature
    # at thermodynamic equilibrium
    ctx0, X0 = model.utils.equilibriate_state(
        context="liquid",
        p=5 * bar,
        T=degC2K(10),
    )
    # attribute this state to all sites
    data.initial_states[...] = X0
    data.initial_contexts[...] = ctx0
