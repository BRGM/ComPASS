import numpy as np
import ComPASS

well = ComPASS.Well()
well.geometry.radius = 0.1
segment = np.reshape(np.arange(10), (-1, 2))
well.geometry.add_segments(segment)
segment = np.reshape(np.arange(120, 140, 2), (-1, 2))
well.geometry.add_segments(segment)
print(well)
print('Well is stopped:', well.stopped)
