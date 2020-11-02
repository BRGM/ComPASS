.. _setting-initial-values:

Setting up initial values
=========================

This script shows several ways to initialize physical values before
running the simulation.

.. literalinclude:: ../test/unit/init_states.py
   :language: python
   :linenos:

You can also reload any previous simulation state from an output directory
using the `simulation.reload_snapshot` method.

.. warning::
    Using the `simulation.reload_snapshot` method
    the mesh and its partition must be exactly the same
    (i.e. you must use exactly the same number of processors for parallel simulations).
    Do NOT forget to reset Dirichlet conditions if necessary.
