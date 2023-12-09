.. meta::
    :scope: version4

.. _setting_linear_solver:

Setting up linear solver options
================================

The way of solving the linear system :math:`Ax=b` is a crucial step of the ComPASS code.
ComPASS can use two types of solvers,
iterative solvers together with a preconditioner,
or direct ones (for small matrix sizes or when iterative solvers do not converge).
By default the linear system is solved using the GMRES procedure
with the :ref:`CPR-AMG <cpramg>` preconditioner.

You can specify the :ref:`linear solver <Linear solvers>` options
with command line options or from the script.

Command line options
--------------------

This uses the `inept paquage <https://pypi.org/project/inept/>`_.

You can try running :download:`doublet_from_options.py <../test/linalg/doublet_from_options.py>`
with different options, for example the following run sequentially with
a Petsc direct solver and display a short view:

.. code:: python

    python3 doublet_from_options.py --lsolver.new.direct True --lsolver.view True

.. literalinclude:: ../test/linalg/doublet_from_options.py
   :language: python
   :linenos:

More informations can be found :ref:`here <Linear solvers>`.

Download file:
:download:`doublet_from_options.py <../test/linalg/doublet_from_options.py>`

Precise the linear solver options in your script
------------------------------------------------

In this case you can execute the file as usual

.. code:: python

    python3 doublet_from_options.py

.. literalinclude:: ../test/linalg/direct_solving.py
   :language: python
   :linenos:

See the :ref:`linear solver section <Linear solvers>` for more details.

Download file:
:download:`direct_solving.py <../test/linalg/direct_solving.py>`

Others examples are located in the `test/linalg` directory.
