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
