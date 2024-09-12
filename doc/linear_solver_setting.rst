.. meta::
    :scope: version5

Setting up linear solver options
================================

The way of solving the linear system :math:`Ax=b` is a
crucial step of the ComPASS code.
ComPASS can use two types of PETSc solvers,
iterative solvers together with a preconditioner,
or direct ones (for small matrix sizes or
when iterative solvers do not converge).
By default the linear system is solved using the GMRES procedure
with the :ref:`CPR-AMG <cpramg>` preconditioner.
If PETSc is not available in the Loaf software brick,
the direct Eigen solver is used instead.

By default the iterative PETSc solver is used in the
:code:`compass-coats.Standard_time_loop` class.
It has default value for the tolerance, the maximum number of iterations...

It is possible:

* to precise to use the direct solver at the creation of the instance of
  the :code:`Standard_time_loop` class

  .. code-block:: python

        time_loop = Standard_time_loop(
            geom=geom,
            model=model,
            scheme=scheme_def,
            data=data,
            output_dir=visu_dir,
            direct_solver=True,
        )

* or to give your own solver at the creation of the instance
  of the :code:`Standard_time_loop` class

  .. code-block:: python

        from loaf.eigen_solver import EigenSolver
        time_loop = Standard_time_loop(
            geom=geom,
            model=model,
            scheme=scheme_def,
            data=data,
            output_dir=visu_dir,
            solver=EigenSolver(),
        )

.. * or to change the value of the iterative solver coefficients,

..   .. code-block:: python

..     time_loop.loop.timestepper.step_solver.linear_solver

.. TODO v5
.. For more details, refer to the :ref:`linear solver section <Linear solvers>`.
