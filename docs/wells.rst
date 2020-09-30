Using wells
===========

Defining wells
--------------

The :func:`simulation.init <ComPASS.simulation.base.init>` has a keyword parameter
that can be used to pass a list of well objects.

A well object is defined in two steps:
    1. its geometry
    2. its role (producer/injector)


The well geometry is defined providing a list of mesh edges in random order.
Edges will be automtically sorted. The only condition is that all edges can be
chained up to the well head that must be unique.
Well edges are defined as a pair of mesh vertices.

Convenience functions are availabe for simple situations:
  - for vertical wells use the :func:`simulation.create_vertical_well <ComPASS.wells.wells.create_vertical_well>`

Once the `well` object is created :
  1. It's a good practice to define a well id which must be a unique integer.
     Though each well will be given a default id, you will be able to reference quickly your well
     with a user defined id :

     .. code-block:: python

        my_well_name = 12345 # <- you choose
        well = simulation.create_vertical_well((0,0))
        well.id = my_well_name

  2. You set the operating conditions (either pressure or flowrate) with :

       .. code-block:: python

          well.operate_on_flowrate = Qw, np.inf # target flow rate in kg/s, threshold pressure

      or:

       .. code-block:: python

          well.operate_on_pressure = pw, Qmax # target pressure Pa, threshold flowrate

      above the threshold pressure/flowrate, the well operating conditions will be switched.

  3. You have to define the role of the well setting it as either a producer with:

     .. code-block:: python

        well.produce()

     or either an injector with:

     .. code-block:: python

        injection_temperature = 300 # K
        well.inject(injection_temperature)


You will find two example in the :ref:`example scripts section<setting_well_transients>`.


Setting well history
--------------------

You can set well transients using the
`simulation.well_production_history` and `simulation.well_injection_history`
functions. You will find two example in the :ref:`example scripts section<setting_well_transients>`.

.. warning::
    For the time being, closed wells are discard during simulation setup.
    `simulation.init` will issue a warning not to distribute closed wells,
    but it's ok to close wells afterwards.


Monitoring well state
---------------------

All well nodes (called perforations) can be acessed and hold the following physical values:
  - pressure at the well node
  - temperature at the well node
  - fluid density at the well node
  - saturations at the well node (an array with number of phases values)
  - pressure drop at the well node
  - molar flowrates at the well node (an array with number of components values)
  - flowing energy at the well nodes

One specific perforation is the well head that can be accessed with
the `simulation.get_wellhead` function, for example:

.. code:: python

    # wid is the well id
    wellhead = simulation.get_wellhead(wid)
    print(f"Well head pressure for well {wid} is: {wellhead.pressure}")


Connections between wells
-------------------------

Connections can be defined between wells so that the well head information
from a given well is made available to another one (whatever the procs that manage the wells).

To connect two wells you give a sequence (list, array...) of pairs `(source, target)`
using the `simulation.add_well_connections` function.

Then the well head information (`molar_flowrate`, `energy_flowrate`, `pressure`, `temperature`)
is made available using the source well id with `simulation.well_connections[source_well_id]`.
For example:

.. code:: python

    wellhead = simulation.well_connections[wid]
    print(f"Well {wid} wellhead pressure is: {wellhead.pressure}")

This can be used to chain well productions using `simulation.standard_loop` iteration callbacks.
The example :download:`chain_random_wells.py <../test/bulk/chain_random_wells.py>` demonstrates
such a use case.

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
how well information can be collected on the master proc and dump at the end of the simulation.

.. note::
    Most of the time you will want to collect well information on the master proc so that
    the simulation script can run both in sequential and parallel.


Error on well settings at convergence
-------------------------------------

The maximum error on well at Newton convergence can be displayed setting:
::

    simulation.newton.check_well_errors_at_convergence = True

Then at the end of each successful Newton loop it will display the maximul on
imposed flowrate and pressure.
