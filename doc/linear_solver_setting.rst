Setting up linear solver options
================================

In the ComPASS code, by default the linear system is solved using the GMRES procedure
with the :ref:`CPR-AMG <cpramg>` preconditioner.

You can modify the linear solver (use a direct one for example)
or specify other linear solver options
with command line (using the `inept paquage <https://pypi.org/project/inept/>`_)
or from the script.

Refer to the :ref:`linear solver section <linear_solvers>` and
to the :ref:`examples section <setting_linear_solver>` for more details.
