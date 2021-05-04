Linear Solvers in ComPASS
=========================

Linear systems in ComPASS are solved using specific objects called
LinearSolver implemented in the ComPASS.linalg package. This document
describes the implementation of linear solvers in ComPASS and sums up a
few strategies for you to manage your own solvers.

Different types of solvers and implementations
----------------------------------------------

In ComPASS there are two kinds of linear solvers, direct and iterative,
each having two different possible implementations. Direct solvers use
an LU decomposition of the matrix to find a solution to the linear
system :math:`Ax=b`. This type of solver can solve any system, but is
not parallelizable and costs about :math:`\frac{2}{3}n^3` operations
(with :math:`n` the number of unknowns/equations), which makes it not
suitable for large problems. It should be used for smaller systems, or
when iterative solvers don’t converge.

Iterative solvers use the GMRES method with preconditioning to
iteratively reduce the distance between an arbitrary first guess and the
solution. The convergence of this type of solvers is highly dependent on
preconditioner quality, meaning that basic preconditioners will lead to
slow convergence or even divergence in case of complex physics, like in
the ComPASS context. The default preconditioning method in ComPASS is
CPR-AMG, an efficient method designed for linear systems in the context
of Darcy flow physics. In general, GMRES with CPR-AMG preconditioning is
the best method in ComPASS simulations. A description of this method can
be found in the docs.

Each of these linear solvers has two different implementations : the
``Legacy`` and the ``Petsc`` versions. The default is still the
``Legacy`` implementation and uses old Fortran routines, while the
latter is a new implementation that uses the petsc4py interface.
Eventually the ``Legacy`` version is bound to disappear, but more tests
have to be made with the newer ``Petsc`` implementation before we can
guarantee accurate and robust results. Advantages and disadvantages of
these two methods are discussed in the next paragraph.

Command line options
--------------------

A default implementation of the new ``PETSc`` version provides the same
the behavior explained in the previous paragraph.

When running a default ``standard_loop(simulation)``, the linear solver
is instanciated internally together with the Newton object that manages
the newton loop, and cannot be accessed from the script that called the
time loop. However, command line options for ComPASS linear solvers are
available for basic usage, and can be used to set the default linear
solver :

-  ``--legacy_linear_solver <legacy/new>`` : Switch between the Fortran
   and the Python implementations (defaults to legacy if not provided)
-  ``--direct_linear_solver`` : Use a direct solver (iterative GMRES algorithm by default)
-  ``--linear_solver_view`` : Display a view of the linear
   solver object at constructor call
-  ``--dump_ls <time>`` : Write the PETSc linear systems (matrix, solution and
   right hand side) in ASCII format files from the end of time step
   ``<time>`` (slow in parallel)
-  ``--dump_ls_binary <time>`` : Write the PETSc linear systems (matrix, solution
   and right hand side) in binary format files from the end of time step
   ``<time>`` (not available in the ``Legacy`` version)
-  ``--abort_on_linear_failure`` : Stop computation if the linear solver fails
-  ``--dump_system_on_linear_failure`` : Dump the linear system if the linear solver fails
-  ``--abort_on_newton_failure`` : Stop computation if the newton loop fails to converge

Examples
~~~~~~~~

Running the following code from file
``test/linalg/solving_from_options.py`` :

.. literalinclude:: ../test/linalg/solving_from_options.py
   :language: python
   :linenos:
   :lines: 3-

with the command :

.. code:: console

   $ python3 solving_from_options.py --linear_solver_version new --direct_linear_solver --linear_solver_view

should display the following, together with the simulation output :

.. code:: console

    LinearSolver object view:
        Direct
        petsc4py new implementation

Instanciate your own linear solver object
-----------------------------------------

The easiest and fastest way of getting your own ``LinearSolver`` object
is to use the factory function from the ``ComPASS.linalg.factory`` file:

.. literalinclude:: ../ComPASS/linalg/factory.py
   :linenos:
   :language: python
   :pyobject: linear_solver
   :end-before: """



Description of the factory function’s arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``simulation`` : The simulation object from which the linear system
   is going to be built
-  ``legacy`` : A boolean switch bewteen the ``Legacy`` and ``Petsc``
   versions
-  ``direct`` : A boolean switch between direct and iterative solver
-  ``activate_cpramg`` : A boolean switch to turn the CPR-AMG
   preconditioner on or off. Defaults to ``None`` in case of direct
   solver, and ``True`` if an iterative solver is chosen.
-  ``cpr_amg_type`` : Select library used for AMG computation in the
   CPR-AMG preconditioner. Possible values : "hypre" (HYPRE BoomerAMG implementation),
   "amg" (PETSc built-in algorithm). Defaults to "hypre".
-  ``tolerance`` : Relative decrease in the residual norm required for
   iterative solver convergence, defaults to 1e-6
-  ``max_iterations`` : Maximum number of iterations accepted before
   iterative solver divergence, defaults to 150
-  ``restart_size`` : Number of iterations at which GMRES restarts,
   defaults to 30
-  ``from_options`` : Enable command line options (set to ``True`` in
   the ``standard_loop()`` call to the factory)

To pass the linear solver object to the simulation loop, you must call
the Newton object constructor. The Newton constructor takes four
arguments : the simulation object, the relative tolerance on residual
decrease to achieve Newton convergence, the maximum number of Newton
iterations, and a linear solver object. The example below shows how to
use a direct solver of the ``Petsc`` kind :

.. literalinclude:: ../test/linalg/direct_solving.py
    :linenos:
    :lines: 3-
    :language: python

Note on the two different versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is deliberately hidden that both implementations actually use the
PETSc library, because we eventually want to get rid of the ``Legacy``
implementation. The ``Legacy`` version uses the Fortran interface, and
the ``Petsc`` version uses the petsc4py bindings. Since both strategies
initialize PETSc, using multiple types of solvers in the same execution
can lead to unexpected behaviour. This means **you should not have
instances of both the** ``Legacy`` **and** ``Petsc`` **versions at the
same time during execution.**

Both version are essentially identical, but there are a few things to
note on the ``Legacy`` version :

-  It has been the default implementation for a while and has been
   tested on many different simulations, so it should be reliable
-  The PETSc objects in the Fortran layer are global, which means that
   it is not possible to have multiple ``Legacy`` linear solver
   instances. A new instanciation will overwrite the settings of the
   first, and both Python instances will refer to the same Fortran Petsc
   objects.

Manipulating linear solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to access an iterative solver’s parameters using the
‘dot’ symbol, and linear solvers also feature a ``__str__`` method, so you can display
information about the object with a ``print()`` statement (or
``mpi.master_print()`` in parallel) :

.. literalinclude:: ../test/linalg/iterative_solving.py
    :linenos:
    :lines: 3-


Output :

::

   LinearSolver object view:
     Iterative
     petsc4py new implementation
     Settings : IterativeSolverSettings(tolerance=1e-07, max_iterations=250, restart_size=50)
     CPR-AMG : activated

Extra docs and examples
-----------------------

A few scripts can be found in the test/linalg directory and can serve as
examples on how to set and use linear solvers in ComPASS. Extra doc on
the CPR-AMG preconditioning method and its implementation in the
``Petsc`` version can be found in the dcumentation.
