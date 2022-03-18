Defining fractures
==================

The fractures are objects of codimension one, they are defined as a set of faces.
The geometry of the fractures is given with the keyword :code:`fracture_faces` in
:func:`simulation.init <ComPASS.simulation.base.init>`.
It can be a list of face global indexes, or a function, or a boolean vector of size :code:`mesh.nb_faces`.

The fracture thickness is up to know a global property which is set
(before calling the :func:`simulation.init <ComPASS.simulation.base.init>` function) as follows:

.. code-block:: python

    fracture_thickness = 1 # m
    simulation.set_fracture_thickness(fracture_thickness)

The permeability, porosity and thermal conductivity are regionalized properties defined
through the
:func:`simulation.init <ComPASS.simulation.base.init>` function and the corresponding keywords:


.. code-block:: python

    simulation.init(
        ...,
        fracture_permeability=k_fracture,
        fracture_porosity=omega_fracture,
        fracture_thermal_conductivity=K_fracture,
        ...,
    )


More information in the :ref:`physical properties section<setting_physical_properties>`
or about the :ref:`physical units<units>`.

Examples with fractures are detailled in the :ref:`example section<classical_doublet_with_frac>`.
