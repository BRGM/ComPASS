.. _physics_section:

Physics
========

Different physics exist, it determines the number of phases and
components and the matrix of presence of the components in the phases.
It also comes with default physical properties (such as the phase densities, viscosities...).

 * linear_water : one phase (by default liquid), one component (by default water).
 * water2ph : two phases (by default liquid and gas), one component (by default water).
 * immiscible2ph : two phases (by default liquid and gas), two components (by default water and air), only water in liquid phase and air in gas phase.
 * diphasic : two phases (by default liquid and gas), two components (by default water and air), all components can be in all phases.


.. code-block:: python

    simulation = ComPASS.load_physics("linear_water")

You can define your own physical properties, refer to :ref:`this section <fluid physical properties>`.
