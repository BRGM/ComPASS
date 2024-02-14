# You can choose here which physics to build and make available in ComPASS
# Available physics are:
#
# * water2ph: subcritical water with 2 phases
# * diphasic: subcritical water with 2 phases and air component, specific
#   atmospheric boundary conditions are made available
# * diphasicCO2: subcritical water with 2 phases and CO2 component
# * linear_water: monophasic water with super simple physical properties
# * immiscible2ph: immiscible liquid water and air in gas phase

option(ComPASS_WITH_water2ph_PHYSICS "build water2ph module" ON)
option(ComPASS_WITH_diphasic_PHYSICS "build diphasic module" ON)
option(ComPASS_WITH_diphasicCO2_PHYSICS "build diphasicCO2 module" ON)
option(ComPASS_WITH_linear_water_PHYSICS "build linear_water module" ON)
option(ComPASS_WITH_immiscible2ph_PHYSICS "build immiscible2ph module" ON)

# Build all physics. If activated, it will override all previous options and
# will include experimental physics.
option(ComPASS_WITH_ALL_PHYSICS "build all physics modules" OFF)
