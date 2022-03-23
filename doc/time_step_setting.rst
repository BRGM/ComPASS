Setting up the time-step options
================================

It is usually useful to manage the time-step such as the initial time-step,
the maximum time-step or the increase (resp. decrease) factor if the time step is
successful (resp. has failed).
The :code:`ComPASS.timestep_management` method gathers these constants as follows:

.. code-block:: python

    from ComPASS.timeloops import TimeStepManager
    ...

    tsmger = TimeStepManager(
        initial_timestep=1 * year,
        minimum_timestep=1e-1, # in s
        maximum_timestep=50.0 * year,
        increase_factor=2.0, # dt^{n+1} = increase_factor * dt^{n} if success
        decrease_factor=0.2, # dt^{n+1} = decrease_factor * dt^{n} if failure
        ...,
    )

This object is then given to the time loop
:func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>` function
with the keyword :code:`time_step_manager`.

.. code-block:: python

    simulation.standard_loop(
        ...,
        time_step_manager=tsmger,
        ...,
    )
