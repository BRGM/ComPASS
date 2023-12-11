.. meta::
    :scope: version5

Numerical scheme
================

Two numerical schemes are currently implemented, VAG and TPFA.

TPFA scheme
-----------

Two point flux approximation scheme: the unknowns are located at the
cells. It is a scheme widly used (for example in Though).

.. code-block:: python

    from compass_coats.schemes import TPFA
    scheme = TPFA()


Careful, the simulation must respect the K-orthogonality.

VAG scheme
----------

This is a numerical scheme where the sites are located at the cell centers
and at the nodes. The cell unknowns can be eliminated using a Schur complement.
Then it ends as a nodal numerical scheme.
It is a scheme very efficient over thetrahedral meshes because they
have much more cells than nodes.

.. code-block:: python

    from compass_coats.schemes import VAG
    scheme = VAG()
