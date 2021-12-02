#pragma once

enum struct Component {
   air = ComPASS_AIR_COMPONENT,
   water = ComPASS_WATER_COMPONENT
};

enum struct Phase { gas = ComPASS_GAS_PHASE, liquid = ComPASS_LIQUID_PHASE };

enum struct Context {
   gas = ComPASS_GAS_CONTEXT,
   liquid = ComPASS_LIQUID_CONTEXT,
   diphasic = ComPASS_DIPHASIC_CONTEXT,
   gas_FF_no_liq_outflow = ComPASS_GAS_FF_NO_LIQ_OUTFLOW_CONTEXT,
   diphasic_FF_no_liq_outflow = ComPASS_DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT,
   diphasic_FF_liq_outflow = ComPASS_DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
};
