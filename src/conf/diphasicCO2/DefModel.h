#pragma once

enum struct Component {
   CO2 = ComPASS_CO2_COMPONENT,
   water = ComPASS_WATER_COMPONENT
};

enum struct Phase { gas = ComPASS_GAS_PHASE, liquid = ComPASS_LIQUID_PHASE };

enum struct Context {
   gas = ComPASS_GAS_CONTEXT,
   liquid = ComPASS_LIQUID_CONTEXT,
   diphasic = ComPASS_DIPHASIC_CONTEXT,
};
