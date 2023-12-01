.. meta::
    :scope: version5

Setting up the time-step manager
================================

It is usually useful to manage the generation of the time-step :math:`dt`.

Two time-step managers are availables, one with a fixed time step
and an adaptative one. If necessary you can implement you own.

Dynamic time-step manager
-------------------------

It is the default time-step manager. There is an *initial_step*; when
the time iteration succeeds, the time step increases by the factor
*increase_factor*; if the time step fails (no convergence of the linear
system solver), it is reduced by the factor *decrease_factor*.

It is possible to adapt the coefficients as follows

.. code-block:: python

    time_loop.loop.timestep_manager.increase_factor = 1.5
    time_loop.loop.timestep_manager.decrease_factor = 0.5
    time_loop.loop.timestep_manager.minimum_step = 100  # seconds

The time steps are also modified to reach the output times and the final time.

Fixed time-step manager
-----------------------

It generates constant time steps, except to reach
the output times and the final time.

Then the fixed time-step manager cannot try an other step if the linear system
solver fails to converge.

In time_loop.run
----------------

You can modify some options in :code:`time_loop.run`:

* *fixed_step*: the time-step manager
  is changed into the :ref:`fixed time-step manager<Fixed time-step manager>`.

* *initial_step*: updates the initial step of the
  :ref:`dynamic time-step manager<Dynamic time-step manager>`.
  It is particularily useful
  when running several times :code:`time_loop.run` to ease the first
  time step convergence.

.. code-block:: python

   solution, tick = time_loop.run(
      fixed_step=0.1 * year,
      final_time=10 * year,
      output_every=10,
   )
