.. meta::
    :scope: version5

Numerical scheme
================

Two numerical schemes are currently implemented, TPFA and VAG.

TPFA scheme
-----------

Two point flux approximation scheme: the sites are located at the
cells and at the Dirichlet faces.
It is a scheme widly used (for example in Though).

.. code-block:: python

    from compass_coats.schemes import TPFA
    scheme_def = TPFA()


Careful, the simulation must respect the K-orthogonality.

VAG scheme
----------

In the Vertex Approximate Gradient scheme the sites are located at the cell centers
and at the nodes. The cell unknowns can be eliminated using a Schur complement.
Then it ends as a nodal numerical scheme.
It is a scheme very efficient over thetrahedral meshes as there is
much more cells than nodes.

.. code-block:: python

    from compass_coats.schemes import VAG
    scheme_def = VAG()

Impact on the script
--------------------

This difference of location of the sites has an impact on the script:

* The :ref:`Dirichlet boundary condition<Dirichlet>` must be set over
  the concerned sites.
  Then with **TPFA**
  the Dirichlet are set over **faces**; whereas with **VAG**
  they are set over **nodes**.

* The :ref:`initial values<Setting up initial values>`
  must be given over all the sites, ie all the **cells** must
  be initialized with **TPFA**; all the **cells and nodes** must be initialized
  with **VAG**.
