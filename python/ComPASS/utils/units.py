#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

km = 1E3 # m

minute = 60 # s
hour = 60 * minute
day = 24 * hour
year = 365.25 * day
ky = 1E3 * year
My = 1E6 * year

bar = 1E5 # Pa
MPa = 1E6 # Pa

ton = 1000 # kg

degC2K = lambda T: T + 273.15
K2degC = lambda T: T - 273.15

time_string = lambda tin: '%10.5g s = %10.5g y' % (tin, tin / year)
