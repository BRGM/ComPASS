import numpy as np
import ComPASS

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
        if well.stopped:
            print('Well is stopped')
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
