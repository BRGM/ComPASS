Load the physics
================

.. ifconfig:: versionlevel <= '4'

  .. include:: physics_summary_v4.rst

.. ifconfig:: versionlevel > '4'

  .. include:: physics_summary_v5.rst

Details about the physics (and their default physical properties)
are in :ref:`this section <Available physics>`.

To load the *diphasic* physics, use:

.. ifconfig:: versionlevel <= '4'

  .. code-block:: python

      simulation = ComPASS.load_physics("diphasic")

  After loading a physics, you can also
  :ref:`set your own physical properties<Fluid properties>`.

.. ifconfig:: versionlevel > '4'

  .. code-block:: python

      from compass_coats.models import Coats
      model = Coats("diphasic")

  You can change the phases, components, contexts names
  to use them in your script and in the visualization.
  **It does not modify the physical properties**, to do so refer to the
  :ref:`the Fluid physical properties section<Fluid properties>`.

  .. code-block:: python

      model = Coats(
        "diphasic",
        # change name of components (only the name, no modif of physical prop!)
        components=["water", "alkane"],
        # change name of phases (only the name, no modif of physical prop!)
        phases=["oil", "liquid"],
        # change name of contexts
        contexts=["oil", "diphasic", "liquid"],
    )
