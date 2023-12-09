.. meta::
    :scope: version4

.. _setting_well_transients:

Setting transient well history
==============================

Well operation history can be turned into time events list
with the `simulation.well_production_history` and `simulation.well_injection_history`
functions.
The events are generated from a sequence (list, array...) of `(time, flowrate)` for producers
and `(time, flowrate, injection temperature)` for injectors.
Then these events are passed to the time loop.

If you have several wells, just concatenate the event lists into a single one.
Events will always be sorted before being processed.

Producer example
----------------

.. literalinclude:: ../test/bulk/well_production_history.py
   :language: python
   :linenos:

Download file:
:download:`well_production_history.py <../test/bulk/well_production_history.py>`

Injector example
----------------

.. literalinclude:: ../test/bulk/well_injection_history.py
   :language: python
   :linenos:

Download file:
:download:`well_injection_history.py <../test/bulk/well_injection_history.py>`
