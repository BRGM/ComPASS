The python "simulation" object
==============================

Most of the parameters used to run a simulation are accessible through
a simulation object. The simulation object is created at the beginning
of the script from one of the available physics, e.g.:
:code:`simulation=ComPASS.load_physics("water2ph")`

.. warning::
    As of today, simulation is a *fake* object that wraps a mix of
    pure python methods and compiled methods that reside in the physics module.

    This is bound, and our mid-term goal is to implement a plain simulation object, with the goal
    to be able to handle several simulation objects at the same time.

    Simulation internals are managed through the :py:mod:`ComPASS._kernel.py` module
    and the :py:mod:`ComPASS.simulation` sub-package.


Accessing the petrophysics properties
-------------------------------------

All petrophysical parameters can be access through a
:code:`simulation.petrophysic()` object that instantiates a wrapper
around several properties. For example to access local properties
you could do:


.. code:: python

    petrophysics = simulation.petrophysics()
    petrophysics.cell_permeability                # a tensor array with shape (nc, 3, 3)
    petrophysics.cell_porosity                    # an array with shape (nc,)
    petrophysics.cell_thermal_conductivity        # a tensor array with shape (nc, 3, 3)
    # And if (and only if) fractures are defined
    petrophysics.fracture_permeability            # an array with shape (nf,)
    petrophysics.fracture_porosity                # an array with shape (nf,)
    petrophysics.fracture_thermal_conductivity    # an array with shape (nf,)


You can change parameters values as the properties are memory
views of the underlying Fortran/C++ arrays.
