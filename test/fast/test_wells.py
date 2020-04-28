#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS


def test_wells():

    # FIXME: this is transitory but we need an eos to acces the well structure
    simulation = ComPASS.load_eos('water2ph')

    W1 = ComPASS.Well()
    W1.geometry.radius = 0.1
    segment = np.reshape(np.arange(10), (-1, 2))
    W1.geometry.add_segments(segment)

    W2 = ComPASS.Well()
    W2.geometry.radius = 0.15
    segment = np.reshape(np.arange(120, 140, 2), (-1, 2))
    W2.geometry.add_segments(segment)

    def display_wells():
        for well in [W1, W2]:
            print()
            print(well)
            if well.closed:
                print('Well is closed')
            else:
                if well.injecting:
                    print('Injection well:')
                elif well.producing:
                    print('Production well:')
                print('    operate on pressure:', well.operate_on_pressure)
                print('    operate on flowrate:', well.operate_on_flowrate)
        print()

    display_wells()

    W1.operate_on_flowrate = 300, 2E5
    W1.produce()
    W2.operate_on_pressure = 4E5, 300
    W2.inject(273.15 + 30)

    display_wells()

    W3 = ComPASS.Well()
    W3.geometry.radius = np.pi
    segments = np.transpose(np.vstack([np.arange(10, 18), np.arange(11, 19)]))
    W3.geometry.add_segments(segments)
    W3.operate_on_pressure = 1E5, 300
    W3.produce()

    print('Well W3:', W3)

    wells = [W1, W2, W3]
    simulation.set_well_geometries(wells)
    simulation.set_well_data(wells)

if __name__== '__main__':
    test_wells()
