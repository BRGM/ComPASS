#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

km = 1e3  # m

minute = 60  # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365.25 * day
ky = 1e3 * year
My = 1e6 * year

bar = 1e5  # Pa
MPa = 1e6  # Pa

ton = 1000  # kg

degC2K = lambda T: T + 273.15
K2degC = lambda T: T - 273.15


def time_string_factory(time_unit, symbol):
    def f(tin):
        return f"{tin:10.5g} s = {tin / time_unit:10.5g} {symbol}"

    return f


time_string = time_string_factory(year, "y")
