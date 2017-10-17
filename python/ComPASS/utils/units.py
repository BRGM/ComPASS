#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#



minute = 60 # s
hour = 60 * minute
day = 24 * hour
year = 365.25 * day

bar = 1E5 # Pa
MPa = 1E6 # Pa

ton = 1000 # kg

degC2K = lambda T: T + 273.15
K2degC = lambda T: T - 273.15
