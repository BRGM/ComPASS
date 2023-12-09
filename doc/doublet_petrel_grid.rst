.. meta::
    :scope: version4

.. _eclipse_grid:

Geothermal doublet on an Eclipse Grid
=====================================

Doublet on a simple cartesian grid
----------------------------------

The grid is extracted using the :code:`PetrelGrid` from the :code:`ComPASS.utils.petrel` module.
The permeability tensors are extracted from the Petrel file using the :code:`PetrelGrid.permeability` property.

.. literalinclude:: ../test/cases/petrel/doublet_on_cartesian_grid_petrel.py
   :language: python
   :linenos:

Download files:
:download:`doublet_on_cartesian_grid_petrel.py <../test/cases/petrel/doublet_on_cartesian_grid_petrel.py>`
:download:`sample.grdecl <../test/cases/petrel/sample.grdecl>`

Doublet on a faulted grid
-------------------------

Where unconformities are present, they are automatically
identified as fracture/faulted surfaces.
Such faces are remeshed relying on CGAL under the hood
(with petrelgridio and pycgal python modules)
to produce conformal meshes.

The number of cells is kept constant, only cell faces are split
to exploit the fact that the VAG scheme is valid on generic
polyhedral meshes.


.. literalinclude:: ../test/cases/petrel/doublet_on_faulted_petrel_grid.py
   :language: python
   :linenos:

Sample petrel grid comes from
`the PyGRDECL project <https://github.com/BinWang0213/PyGRDECL/tree/master/ExampleData>`_.

Download files:
:download:`doublet_on_faulted_petrel_grid.py <../test/cases/petrel/doublet_on_faulted_petrel_grid.py>`
:download:`Simple20x20x5_Fault.grdecl <../test/cases/petrel/Simple20x20x5_Fault.grdecl>`
