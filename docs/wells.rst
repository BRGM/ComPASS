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

This can be used to chain well productions using `simulation.standard_loop` iteration callbacks.
The example :download:`chain_random_wells.py <../test/bulk/chain_random_wells.py>` demonstrate
such a use case.

.. warning::
    Doing so wells are chained but not coupled. So the simulation result will strongly depend
    on the timestep (do not take too big a timestep).
