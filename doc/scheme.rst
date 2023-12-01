.. meta::
    :scope: version5

Numerical scheme
================

Two numerical schemes are currently implemented, VAG and TPFA.

TPFA
----

The simulation can use the two point flux approximation scheme:

.. code-block:: python

    from compass_coats.schemes import TPFA
    scheme = TPFA()


Careful, the simulation must respect the K-anisotropy.

VAG
----

This is a numerical scheme where the sites are located at the cell centers
and at the nodes. The cell unknowns can be eliminated using a Schur complement.
Then it ends as a nodes centered numerical scheme.

.. code-block:: python

    from compass_coats.schemes import VAG
    scheme = VAG()
