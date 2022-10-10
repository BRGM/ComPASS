.. _bc:

Setting up the boundary conditions
==================================

Dirichlet nodes
---------------

The set of Dirichlet nodes can be identified in the
:func:`simulation.init<ComPASS.simulation.base.init>`
method with the keyword :code:`set_dirichlet_nodes`, for example:

.. code-block:: python

    simulation.init(
        mesh = grid,
        ...,
        set_dirichlet_nodes=simulation.vertical_boundaries(grid),
        ...,
    )

The function :func:`simulation.reset_dirichlet_nodes` modifies the
set of Dirichlet nodes.
It is useful to modify the set of Dirichlet nodes, for example
between two calls of the time loop :func:`simulation.standard_loop`
function.
:func:`simulation.reset_dirichlet_nodes` is also convenient to identify
the set of Dirichlet nodes after calling the
:func:`simulation.init<ComPASS.simulation.base.init>` function, which has the
advantage to deal with local arrays (after the partition).

.. code-block:: python

    simulation.init(
        mesh = grid,
        ...,
    )
    simulation.reset_dirichlet_nodes(on_vertical_boundaries(grid))

It is also possible to distinguish between the Temperature and the Pressure
set of Dirichlet nodes, for example:

.. code-block:: python

    simulation.init(
        mesh = grid,
        ...,
        set_pressure_dirichlet_nodes=simulation.vertical_boundaries(grid),
        set_temperature_dirichlet_nodes=simulation.all_boundaries(grid),
        ...,
    )

The Dirichlet nodes must be initialized using the set of states contained in
:func:`simulation.dirichlet_node_states`. Different ways to initialize the
nodes are presented in :ref:`this script example<setting-initial-values>`,
follows an example:

.. code-block:: python

    # build a liquid state at fixed pressure and temperature at equilibrium
    X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
    # attribute this state to all the Dirichlet nodes
    simulation.dirichlet_node_states().set(X0)

.. warning::
    Must be at thermodynamic equilibrium, the Dirichlet nodes remain in the
    computation of the residual !


.. _neumann_faces_bc:

Neumann faces
-------------

Neumann boundary condition are set via the
:func:`simulation.set_Neumann_faces` function
using a special ComPASS object called
:code:`ComPASS.NeumannBC` (containing the heat flux in
:math:`W.m^{-2}` and/or the molar flux in
:math:`mol.m^{-2}.s^{-1}`).
Careful: then molar flux must be defined for each component !

**Remark**: if using the *water2ph* physics, the :ref:`system writes the
mass balance equation instead of the molar balance equation<water2ph_equations>`,
then :code:`ComPASS.NeumannBC().molar_flux` must be initialized
with the mass flux expressed in :math:`kg.m^{-2}.s^{-1}`.

.. code-block:: python

    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux # in W/m^2 = J/m^2/s
    Neumann.molar_flux[:] = Qm # one value by component (in mol/m^2/s)
    face_centers = simulation.face_centers()
    simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)

Or in a synthetic way:

.. code-block:: python

    Neumann = ComPASS.NeumannBC(Qm, bottom_heat_flux)
    face_centers = simulation.face_centers()
    simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)

.. _frac_edges_bc:

Neumann fracture edges
----------------------

It is also common to set Neumann boundary condition only on some fracture
edges, the function :func:`simulation.set_Neumann_fracture_edges` set the
values and the edges. The function :func:`simulation.find_fracture_edges` is
useful in this case to extract a set of fracture edges from a set
of faces (or a mask). Follows an example:

.. code-block:: python

    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = Qm # one value by component
    Neumann.heat_flux = bottom_heat_flux
    face_centers = simulation.face_centers()
    where = (
        (np.abs(face_centers[:, 0]) < 0.25 * dx)
        & (np.abs(face_centers[:, 1]) <= 0.2 * Ly)
        & (np.abs(face_centers[:, 2]) <= 0.25 * dz)
    )
    left_fracture_edges = simulation.find_fracture_edges(where)
    simulation.set_Neumann_fracture_edges(left_fracture_edges, Neumann)

.. _atmBC:

Atmospheric boundary condition
------------------------------

An atmospheric boundary condition has been developed and implemented in the
ComPASS code. It is implemented only with the *diphasic* physic.
The first step is to identify the set of boundary faces using
the function :func:`simulation.set_freeflow_faces`:

.. code-block:: python

    simulation = ComPASS.load_physics("diphasic")
    ...
    fc = simulation.compute_face_centers()
    simulation.set_freeflow_faces(on_zmax(grid)(fc))

Then it is necessary to initialize the boundary nodes, to do so it is possible
to retrieve the set of nodes using :func:`simulation.get_freeflow_nodes`.

Careful: the set of context are distinct !
They are *gas_FF_no_liq_outflow*, *diphasic_FF_no_liq_outflow* and
*diphasic_FF_liq_outflow*.

:ref:`This script example<setting-initial-values>` presents different ways to
initialize the nodes. Follows a synthetic way using the :code:`set` function:

.. code-block:: python

    is_ff = simulation.get_freeflow_nodes()  # array of bool of size n_vertices
    X_top = simulation.build_state(
        simulation.Context.gas_FF_no_liq_outflow, p=Patm, T=Tinit, Cag=0.99,
    )
    simulation.node_states().set(is_ff, X_top)

.. _far_field_atmBC:

Change the far-field values
...........................

To modify the far-field values, the :func:`simulation.set_atm_...` are useful
if the values are constant in space, such as:

.. code-block:: python

    simulation.set_atm_temperature(Tatm)
    simulation.set_atm_temperature(Patm)
    simulation.set_atm_rain_flux(-3.2e-2)  # mol/m^2/s : 0 by default

If the far-field values are not constant in space, the object
:code:`simulation.freeflow_node_states` can be accessed and modified:

.. code-block:: python

    gasPhase = simulation.phase_index(simulation.Phase.gas)
    liquidPhase = simulation.phase_index(simulation.Phase.liquid)

    is_ff = simulation.get_freeflow_nodes() # array of bool of size n_vertices
    ff_ns = simulation.freeflow_node_states() # Far-field values arrays
    vertices = simulation.vertices()
    zmid = (z_max - z_min) / 2.0
    top_vertices = vertices[:,2] > zmid
    ff_ns.p[is_ff] = patm # scalar value : only gas pressure
    ff_ns.T[is_ff,gasPhase] = Tatm
    ff_ns.T[is_ff,liquidPhase] = Train
    # the following molar flux contains one value for each component (air=0, water=1)
    ff_ns.imposed_flux[np.logical_and(is_ff, top_vertices), waterComponent] = -3.0e-2
    ff_ns.Hm[is_ff, :] = 0.0
    ff_ns.HT[is_ff] = 0.0

.. warning::

    The far-field values are distinct from the boundary nodes values !
    The far-field values are accessed with
    :code:`simulation.freeflow_node_states`
    whereas the boundary values are accessed with the usual object
    :code:`simulation.node_states`.
