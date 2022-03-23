.. _setting_newton:

Setting up the Newton algorithm options
=======================================

It is often necessary to adapt the Newton algorithm parameters to obtain the
convergence of the Newton algorithm. The arguments must be the simulation object,
the tolerance applied to consider the convergence, the maximum number of iterations
and a linear solver object (:code:`simulation.linear_solver()` provides the default one).
Refer to :ref:`this section<linear_solvers_script>` to modify the linear solver options.
If the maximum number of iterations is achieved, the Newton algorithm has failed and
the time loop may decide to try with a smaller time step.

.. code-block:: python

    from ComPASS.newton import Newton
    ....
    newton = Newton(simulation, 1e-6, 20, simulation.linear_solver())

Then give it to the time loop as follows:

.. code-block:: python

    simulation.standard_loop(
        ...,
        newton=newton,
        ...,
    )

.. warning::

    The linear solver tolerance must be smaller than the Newton tolerance.
