.. _wells_introduction:

Using wells
===========

Defining wells
--------------

The :func:`simulation.init <ComPASS.simulation.base.init>` function has a keyword parameter
(*wells*) that is used to pass a function which creates a list of well objects.

A well object is created using its geometry providing a list of
oriented edges (from wellhead downwards) in random order, with
the :func:`simulation.create_well_from_segments <ComPASS.wells.wells.create_well_from_segments>`
function. Edges will be automatically sorted. The only condition is that all edges can be
chained up to the well head that must be unique.
Well edges are defined as a pair of mesh vertices oriented from wellhead downwards.

Convenience functions are availabe for simple geometries
  - for vertical wells, create the well geometry using the
    :func:`well = simulation.create_vertical_well((Wx,Wy)) <ComPASS.wells.wells.create_vertical_well>` function,
  - for a well described by a list of ordered nodes (describing the well from top to bottom), create the well geometry using the
    :func:`well = simulation.create_single_branch_well(nodes) <ComPASS.wells.wells.create_single_branch_well>` function.

Once the `well` object is created :
  1. It's a good practice to define a well id which must be a unique integer.
     Though each well will be given a default id, you will be able to reference quickly your well
     with a user defined id :

     .. code-block:: python

        my_well_name = 12345 # <- you choose
        well = simulation.create_vertical_well((0,0)) # define the geometry
        well.id = my_well_name

  2. You set the operating conditions (either pressure or flowrate) with :

     .. code-block:: python

          well.operate_on_flowrate = Qw, np.inf # target flow rate in kg/s, threshold pressure

     or:

     .. code-block:: python

          well.operate_on_pressure = pw, Qmax # target pressure Pa, threshold flowrate

     Above the threshold pressure/flowrate, the well operating conditions will be switched.

  3. You have to define the role of the well, setting it either as a producer with:

     .. code-block:: python

        well.produce()

     either as an injector with:

     .. code-block:: python

        injection_temperature = 300 # K
        well.inject(injection_temperature)


You will find two examples in the :ref:`example scripts section<classical_doublet>`.


Setting well history
--------------------

You can set well transients using the
`simulation.well_production_history` and `simulation.well_injection_history`
functions. You will find two examples in the :ref:`example scripts section<setting_well_transients>`.

.. warning::
    For the time being, closed wells are discarded during simulation setup.
    `simulation.init` will issue a warning not to distribute closed wells,
    but it's possible to close a well after `simulation.init` to start the simulation
    with a closed well.


Monitoring well state
---------------------

All well nodes (called perforations) can be acessed and hold the following physical values:
  - pressure at the well node
  - temperature at the well node
  - fluid density at the well node
  - saturations at the well node (an array with number of phases values)
  - pressure drop at the well node
  - molar flowrates at the well node (an array with number of components values)
  - flowing energy at the well node

One specific perforation is the well head that can be accessed with
the `simulation.get_wellhead` function, for example:

.. code:: python

    # wid is the well id
    wellhead = simulation.get_wellhead(wid)
    print(f"Well head pressure for well {wid} is: {wellhead.pressure}")

To access all perforations state you can use
the `simulation.get_well_perforations_state` function.
Then, there is no array wrapper to access underlying
property yet (this is a work in progress cf. issue
`298 <https://gitlab.inria.fr/charms/ComPASS/-/issues/298>`_
). But you can easily build a copy :

.. code:: python

    # wid is the well id
    perfs = simulation.get_well_perforations_state(wid)
    p = np.array([perf.pressure for perf in perfs])


Connections between wells
-------------------------

Connections can be defined between wells so that the well head information
from a given well is made available to another one (whatever the procs that manage the wells
when running in parallel).

To connect two wells you give a sequence (list, array...) of pairs `(source, target)`
using the
:func:`simulation.add_well_connections <ComPASS.wells.connections.add_well_connections>` function.

Then the well head information (`molar_flowrate`, `energy_flowrate`, `pressure`, `temperature`)
is made available using the source well id with `simulation.well_connections[source_well_id]`.
For example:

.. code:: python

    wellhead = simulation.well_connections[wid]
    print(f"Well {wid} wellhead pressure is: {wellhead.pressure}")

This can be used to chain well productions using *iteration_callbacks*
in the :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>` function.
The example :download:`chain_random_wells.py <../test/bulk/chain_random_wells.py>` demonstrates
such a use case, where the flowrate and the temperature in each injector well depend on the
connected productor well.

.. warning::
    Doing so wells are chained but not coupled. So the simulation result will strongly depend
    on the timestep (do not take too big a timestep).

You can also make some wells available on a specific processor using the `proc_requests`
keyword of the `simulation.add_well_connections` function. For example:

.. code:: python

    simulation.add_well_connections(proc_requests=[
        (0, [1, 2, 6]), # will make wells 1, 2 and 6 available on proc 0
        (1, [0, 1]), # will make wells 0 and 1 available on proc 1
    ])

The example :download:`chain_random_wells.py <../test/bulk/chain_random_wells.py>` also demonstrates
how well information can be collected on the master proc and dumped at the end of the simulation.

.. note::
    Most of the time you will want to collect well information on the master proc so that
    the simulation script can run both in sequential and parallel.


Error at Newton convergence on well
-----------------------------------

The maximum error on well at Newton convergence can be displayed setting:
::

    simulation.newton.check_well_errors_at_convergence = True

Then at the end of each successful Newton loop it will display the maximum error on
imposed flowrate and pressure.
