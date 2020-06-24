#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# This comes from:
# @techreport{Thiery2012,
# author = {Thi{\'{e}}ry, Dominique},
# file = {:D$\backslash$:/work/data/biblio/Thi{\'{e}}ry - 2012 - Validation des calculs de transport d'{\'{e}}nergie dans le logiciel MARTHE.pdf:pdf},
# institution = {BRGM},
# keywords = {charms cas test},
# mendeley-tags = {charms cas test},
# title = {{Validation des calculs de transport d'{\'{e}}nergie dans le logiciel MARTHE}},
# url = {http://www.brgm.fr/sites/default/brgm/logiciels/marthe/nt{\_}eau{\_}2011{\_}02{\_}marthe{\_}valid{\_}thermique.pdf},
# year = {2012}
# }

import numpy as np

# from ComPASS.utils.units import deg2C
# from ComPASS.timeloops import standard_loop, TimeStepManager

# degK2C = lambda T: T - 273.15
# degC2K = lambda T: T + 273.15


class MA2_analytical:
    def __init__(
        self,
        omega=1.0,  # porosity (fluid only)
        b=1.0,  # wall thickness (m)
        L=4.0,  # wall length (m)
        QE=30.0,  # lateral heat flux (W/m2)
        TL=25.0,  # initial and imposed right temperature
        specific_heat=1.8e6,  # specific heat (J/m3/K)
        thermal_conductivity=5.5,  # thermal conductivity (W/m/K)
    ):
        self.omega = omega
        self.b = b
        self.L = L
        self.QE = QE
        self.TL = TL
        self.specific_heat = specific_heat
        self.thermal_conductivity = thermal_conductivity
        self.D = thermal_conductivity / specific_heat

    def __call__(self, x, t, precision=1e-10):
        assert precision > 0
        D = self.D
        L = self.L
        TL = self.TL
        QE = self.QE
        thcond = self.thermal_conductivity
        td = np.reshape((D * np.asarray(t, dtype=np.double)) / (4 * L ** 2), (-1, 1))
        xd = np.ravel(np.asarray(x, dtype=np.double)) / L
        Td = 1 - np.tile(xd, (td.shape[0], 1))
        # Solution by D. Thiéry - there are errors here
        n = int((1 + np.sqrt(8 / (np.pi ** 2 * precision))) / 2) + 1
        for k in range(n):
            w = (2 * (k + 1) - 1) * np.pi
            Td += (8 / w ** 2) * np.cos(0.5 * w * xd) * np.exp(-(w ** 2) * td)

        return TL + ((QE * L) / thcond) * Td
