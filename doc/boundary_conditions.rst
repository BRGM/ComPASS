.. meta::
    :scope: version5

Setting up boundary conditions
==============================

The boundary condition is set in the
:code:`data.boundary_conditions` object.
The location of the boundary sites relies on the initialized values.
Moreover, the boundary condition objects depend on the numerical scheme,
in :ref:`TPFA<TPFA scheme>` it concerns the boundary faces;
in :ref:`VAG<VAG scheme>` the boundary nodes.
Thus the boundary sites are the intersection between the possible boundary
sites of the numerical scheme and the initialized values.

Thus, it is possible to initialize more objects than the necessary ones
to have a generic script for various numerical schemes.
Utilities for cartesian meshes are in
`geom_traits <https://gitlab.com/compass/compass-v5/geom-traits/-/blob/main/src/geom_traits/grid.py?ref_type=heads>`_
to get the boundary objects.

Dirichlet
---------

It is necessary to set both
:code:`data.boundary_conditions.Dirichlet_states` and
:code:`data.boundary_conditions.Dirichlet_contexts`.
The context is directly related to the Coats formulation, it defines the
present phases at a specific time step and a specific sites.


.. code-block:: python

    data = model.new_data(geom.mesh.dimension)
    # build a liquid state at fixed pressure and temperature
    # at thermodynamic equilibrium
    top_ctx, Xtop = model.utils.equilibriate_state(
        context="liquid",
        p=20 * bar,
        T=degC2K(50),
    )
    # top contains nodes, edges, faces
    # thus possible to use with all numerical scheme
    top = geom_traits.top_boundary(geom)
    # attribute this state to all the top objects
    data.boundary_conditions.Dirichlet_states[top] = Xtop
    data.boundary_conditions.Dirichlet_contexts[top] = top_ctx

.. warning::

    Must be at thermodynamic equilibrium, the Dirichlet sites remain in the
    computation of the residual !


Neumann flux
------------

The Neumann fluxes are set via the
:code:`data.boundary_conditions.Neumann_flux` object.
They are always imposed at faces (no disctinction between numerical schemes).
**The Neumann flux is a surface flux, it will be multiplied by faces measures.**

Careful: **the molar flux must be defined for each component !**

**Remark**: if using the *water2ph* physics, the :ref:`system writes the
mass balance equation instead of the molar balance equation<water2ph_equations>`,
then :code:`Neumann_flux.molar_flux` must be initialized
with the mass flux expressed in :math:`kg.m^{-2}.s^{-1}`.

.. code-block:: python

    model = Coats("diphasic")
    data = model.new_data(geom.mesh.dimension)
    Neumann_flux = data.boundary_conditions.Neumann_flux.dtype()
    Neumann_flux.energy = bottom_heat_flux  # in W/m^2 = J/m^2/s
    Neumann_flux.molar[model.components["air"]] = 0  # in mol/m^2/s
    Neumann_flux.molar[model.components["water"]] = qmol  # in mol/m^2/s
    bottom_faces = geom_traits.bottom_boundary_faces(geom)
    # init all the bottom faces with the same Neumann flux
    data.boundary_conditions.Neumann_flux[bottom_faces] = Neumann_flux
