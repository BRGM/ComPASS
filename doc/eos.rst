Equations of state
==================

Different *equations of state* exist. It determines the possible phases and
components and the matrix of presence of the components in the phases.
It also fix some physical laws (phase density, ...).
 - linear_water : one phase (liquid), one component (water).
 - water2ph : two phases (liquid and gas), one component (water).
 - immiscible2ph : two phases (liquid and gas), two components (water and air), only water in liquid phase and air in gas phase.
 - diphasic : two phases (liquid and gas), two components (water and air), all components can be in all phases.

.. code-block:: python

    simulation = ComPASS.load_eos("linear_water")
