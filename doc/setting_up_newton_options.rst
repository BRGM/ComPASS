.. meta::
    :scope: version5

Setting up the Newton options
=============================

It is often necessary to adapt the Newton algorithm parameters to obtain the
convergence of the Newton algorithm. Those parameters are
the tolerance applied to consider the convergence, and
the maximum number of iterations.

.. code-block:: python

    # Newton options
    time_loop.loop.timestepper.step_solver.tolerance = 1e-6
    time_loop.loop.timestepper.step_solver.maxiter = 25

If the maximum number of iterations is reached, the Newton algorithm
has failed and the time loop may continue with a smaller time step
depending on the :ref:`time-step manager <Setting up the time-step manager>`.

.. warning::

    The linear solver tolerance must be smaller than the Newton tolerance.
