.. meta::
    :scope: version4

Setting up initial inputs
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
    The default is that Dirichlet node states are updated according to the reloaded snapshot.
    You can set the `reset_dirichlet` argument to `False` if you do not want this behavior.

Download file:
:download:`init_states.py <../test/unit/init_states.py>`
