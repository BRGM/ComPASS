import numpy as np
import ComPASS

well = ComPASS.Well()
well.geometry.radius = 0.1
well.geometry.add_branch(np.arange(10))
well.geometry.add_branch(np.arange(120, 130, 2))

print(well)
print('Well is stopped:', well.stopped)
