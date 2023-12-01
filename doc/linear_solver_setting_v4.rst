.. meta::
    :scope: version4

Setting up linear solver options
================================

The way of solving the linear system :math:`Ax=b` is a crucial step of the ComPASS code.
ComPASS can use two types of solvers,
iterative solvers together with a preconditioner,
or direct ones (for small matrix sizes or when iterative solvers do not converge).
By default the linear system is solved using the GMRES procedure
with the :ref:`CPR-AMG <cpramg>` preconditioner.

You can specify the :ref:`linear solver <Linear solvers>` options
with :ref:`command line options <command_line_options>` or from the script.


Precise the linear solver options in your script
------------------------------------------------

In this case, create a `linear_solver` object with the option `direct=True`
and give it to the Newton object. The linear solver is an option of the Newtown loop.
Then precise the Newton in the :func:`simulation.standard_loop <ComPASS.timeloops.standard_loop>`.

.. code-block:: python

    from ComPASS.linalg.factory import linear_solver
    from ComPASS.newton import Newton

    lsolver = linear_solver(simulation, direct=True)
    newton = Newton(simulation, 1e-5, 8, lsolver)

    simulation.standard_loop(
        ...
        newton=newton,
        ...
    )


Command line options
--------------------

This uses the `inept paquage <https://pypi.org/project/inept/>`_.

Precise in your script that the linear solver will be set from command line options:

.. code-block:: python

    from ComPASS.linalg.factory import inept_linear_solver

    command_line_lsolver = inept_linear_solver(simulation)
    newton = Newton(simulation, 1e-5, 8, command_line_lsolver)

    simulation.standard_loop(
        ...
        newton=newton,
        ...
    )


You can define :ref:`various options <command_line_options>` when running your script,
for example the following runs with
a Petsc direct linear solver and displays a short view:

.. code:: python

    python3 my_script.py --lsolver.new.direct True --lsolver.view True


For more details, refer to the :ref:`linear solver section <Linear solvers>`.
