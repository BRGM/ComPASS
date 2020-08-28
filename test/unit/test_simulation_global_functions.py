#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS

simulation = ComPASS.load_eos("water2ph")

print(f"gravity:", simulation.get_gravity())
simulation.set_gravity(10)
assert simulation.get_gravity() == 10

print(f"fracture thickness:", simulation.get_fracture_thickness())
simulation.set_fracture_thickness(0.1)
assert simulation.get_fracture_thickness() == 0.1
