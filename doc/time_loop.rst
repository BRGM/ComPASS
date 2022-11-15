.. _time_loop_exe:

Time loop execution
===================

The time loop is usualy executed using the
:func:`standard_loop<ComPASS.simulation.timeloops.standard_loop>`
function.

It allows lots of optional values to set output options or
time-step specifications or callbacks...

The minimal use is to specify:

* the *final_time* or the maximal number of iterations *nitermax*,

* the *initial_timestep* or the value of the *fixed_timestep*

.. code-block:: python

   simulation.standard_loop(
      initial_timestep=100 * day,
      final_time=1 * year,
   )
