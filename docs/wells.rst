Wells
=====

The maximum error on well at Newton convergence can be displayed setting:
::

    simulation.newton.check_well_errors_at_convergence = True

Then at the end of each successful Newton loop it will display the maximul on
imposed flowrate and pressure.


You can set well transients using the
`simulation.well_production_history` and `simulation.well_injection_history`
functions. You will find two example in the :ref:`example scripts section<setting_well_transients>`.

.. warning::
    For the time being, closed wells are discard during simulation setup.
    `simulation.init` will issue a warning not to distribute closed wells,
    but it's ok to close wells afterwards.


Connections between wells
-------------------------

Connections can be defined between wells so that the well head information
from a given well is made available to another one (whatever the procs that manage the wells).

To connect two wells you give a sequence (list, array...) of pairs `(source, target)`
using the `simulation.add_well_connections` function.

Then the well head information (`mass_flowrate`, `energy_flowrate`, `pressure`, `temperature`)
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
