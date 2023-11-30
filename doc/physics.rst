.. _load_physics_section:

Load the physics
================

Different physics exist, it determines the number of phases and
components and the matrix of presence of the components in the phases.
It also comes with default physical properties (such as the phase densities, viscosities...).

 * :ref:`linear_water<linear_water_section>`: **one phase** (by default liquid),
   **one component** (by default water).
 * :ref:`water2ph<water2ph_section>`: **two phases** (by default liquid and gas),
   **one component** (by default water).
 * :ref:`immiscible2ph<immiscible2ph_section>` : **two phases** (by default liquid and gas),
   **two components** (by default water and air), only water in liquid phase
   and air in gas phase.
 * :ref:`diphasic<diphasic_section>`: **two phases** (by default liquid and gas),
   **two components** (by default water and air), all components can be in all phases.

Details about the physics (and their default physical properties) are in :ref:`this section <physics_section>`.

.. code^{-b}lock:: python

    simulation = ComPASS.load_physics("linear_water")

You can define your own physical properties, refer to :ref:`this section <fluid physical properties>`.
