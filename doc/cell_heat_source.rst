.. meta::
    :scope: version5

Setting up cell heat sources
============================

Through the :code:`data.cell_heat_source` object, it is possible to set
a heat source over any set of cells. It is a volumic value,
it is latter multiplied by the cell volume.

.. code-block:: python

    # define the cell domain
    xyz = Data(geom.cell_centers)
    heating_cells = (xyz[0] < L / 3) & (xyz[1] < L / 2) & (xyz[2] < L / 4)
    # set the heat source
    data.cell_heat_source[heating_cells] = 3.0  # W/m^3 = J/m^3/s

An exemple script is in `test_heat_source.py. <https://gitlab.com/compass/compass-v5/compass-coats/-/blob/main/test/test_heat_source.py?ref_type=heads>`_
