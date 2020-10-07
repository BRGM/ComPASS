#include <pybind11/numpy.h>

#include "Model_wrappers.h"
#include "StateObjects.h"

static_assert(X::Model::np == 2, "Wrong numpber of phases.");
static_assert(X::Model::nc == 2, "Wrong numpber of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS == 6, "Wrong number of contexts.");

enum struct Component {
   air = ComPASS_AIR_COMPONENT,
   water = ComPASS_WATER_COMPONENT
};

enum struct Phase { gas = ComPASS_GAS_PHASE, liquid = ComPASS_LIQUID_PHASE };

// FIXME: assuming liquid phase is the latest phase
// FIXME: cpp starts at 0 where fortran at 1
static_assert(static_cast<int>(Phase::gas) == 1, "Wrong gas phase id.");
static_assert(static_cast<int>(Phase::liquid) == 2, "Wrong liquid phase id.");

enum struct Context {
   gas = ComPASS_GAS_CONTEXT,
   liquid = ComPASS_LIQUID_CONTEXT,
   diphasic = ComPASS_DIPHASIC_CONTEXT,
   gas_FF_no_liq_outflow = ComPASS_GAS_FF_NO_LIQ_OUTFLOW_CONTEXT,
   diphasic_FF_no_liq_outflow = ComPASS_DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT,
   diphasic_FF_liq_outflow = ComPASS_DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
};

// Fortran functions
extern "C" {
void FluidThermodynamics_molar_density(int, double, double, const double *,
                                       const double *, double &, double &,
                                       double &, double *, double *);
void FluidThermodynamics_molar_enthalpy(int, double, double, const double *,
                                        const double *, double &, double &,
                                        double &, double *, double *);
void FluidThermodynamics_dynamic_viscosity(int, double, double, const double *,
                                           const double *, double &, double &,
                                           double &, double *, double *);
void FluidThermodynamics_Psat(double, double &, double &);
void FluidThermodynamics_Tsat(double, double &, double &);
}

void init_model() {}

void finalize_model() {}

template <int PHASE>
inline auto saturate_phase() {
   auto S = std::array<double, X::Model::np>{};  // zero initialization
   static_assert(PHASE > 0, "Inconsistent phase id.");
   static_assert(PHASE <= X::Model::np, "Inconsistent phase id.");
   S[PHASE - 1] = 1;  // FIXME: cpp starts at 0 where fortran at 1
   assert(std::accumulate(begin(S), end(S), 0) == 1);
   return S;
}

template <int PHASE>
inline double phase_molar_density(double p, double T,
                                  py::array_t<double, py::array::c_style> &C) {
   double xsi, dxsidp, dxsidT;
   auto S = saturate_phase<PHASE>();
   double dxsidC[X::Model::nc] = {0};
   auto dxsidS = std::array<double, X::Model::np>{};  // zero initialization
   FluidThermodynamics_molar_density(PHASE, p, T, C.data(), S.data(), xsi,
                                     dxsidp, dxsidT, dxsidC, dxsidS.data());
   return xsi;
}

inline double liquid_molar_density(double p, double T,
                                   py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density<static_cast<int>(Phase::liquid)>(p, T, C);
}

inline double gas_molar_density(double p, double T,
                                py::array_t<double, py::array::c_style> &C) {
   return phase_molar_density<static_cast<int>(Phase::gas)>(p, T, C);
}

template <int PHASE>
inline double phase_molar_enthalpy(double p, double T,
                                   py::array_t<double, py::array::c_style> &C) {
   double h, dhdp, dhdT;
   auto S = saturate_phase<PHASE>();
   double dhdC[X::Model::nc] = {0};
   auto dhdS = std::array<double, X::Model::np>{};  // zero initialization
   FluidThermodynamics_molar_enthalpy(PHASE, p, T, C.data(), S.data(), h, dhdp,
                                      dhdT, dhdC, dhdS.data());
   return h;
}

inline double liquid_molar_enthalpy(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy<static_cast<int>(Phase::liquid)>(p, T, C);
}

inline double gas_molar_enthalpy(double p, double T,
                                 py::array_t<double, py::array::c_style> &C) {
   return phase_molar_enthalpy<static_cast<int>(Phase::gas)>(p, T, C);
}

template <int PHASE>
inline double phase_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   double mu, dmudp, dmudT;
   auto S = saturate_phase<PHASE>();
   double dmudC[X::Model::nc] = {0};
   auto dmudS = std::array<double, X::Model::np>{};  // zero initialization
   FluidThermodynamics_dynamic_viscosity(PHASE, p, T, C.data(), S.data(), mu,
                                         dmudp, dmudT, dmudC, dmudS.data());
   return mu;
}

inline double gas_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity<static_cast<int>(Phase::gas)>(p, T, C);
}

inline double liquid_dynamic_viscosity(
    double p, double T, py::array_t<double, py::array::c_style> &C) {
   return phase_dynamic_viscosity<static_cast<int>(Phase::liquid)>(p, T, C);
}

inline double Psat(double T) {
   double result;
   double dPsatdT;
   FluidThermodynamics_Psat(T, result, dPsatdT);
   return result;
}

inline double Tsat(double p) {
   double result;
   double dTsatdp;
   FluidThermodynamics_Tsat(p, result, dTsatdp);
   return result;
}

void add_specific_model_wrappers(py::module &module) {
   module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
   module.def("gas_molar_enthalpy", py::vectorize(gas_molar_enthalpy));
   module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
   module.def("gas_molar_density", py::vectorize(gas_molar_density));
   module.def("liquid_dynamic_viscosity",
              py::vectorize(liquid_dynamic_viscosity));
   module.def("gas_dynamic_viscosity", py::vectorize(gas_dynamic_viscosity));
   module.def("Psat", py::vectorize(Psat));
   module.def("Tsat", py::vectorize(Tsat));

   py::enum_<Component>(module, "Component")
       .value("air", Component::air)
       .value("water", Component::water);

   py::enum_<Context>(module, "Context")
       .value("gas", Context::gas)
       .value("liquid", Context::liquid)
       .value("diphasic", Context::diphasic)
       .value("gas_FF_no_liq_outflow", Context::gas_FF_no_liq_outflow)
       .value("diphasic_FF_no_liq_outflow", Context::diphasic_FF_no_liq_outflow)
       .value("diphasic_FF_liq_outflow", Context::diphasic_FF_liq_outflow);

   py::enum_<Phase>(module, "Phase")
       .value("gas", Phase::gas)
       .value("liquid", Phase::liquid);

   module.def(
       "build_state",
       [](Context context, double p, double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> S,
          py::array_t<double, py::array::c_style | py::array::forcecast> C,
          py::object outflow_mass_flowrates) {
          if (!(S.ndim() == 1 && S.size() == X::Model::np))
             throw std::runtime_error(
                 "Saturations vector has wrong dimensions.");
          if (!(C.ndim() == 2 && C.shape(0) == X::Model::np &&
                C.shape(1) == X::Model::nc))
             throw std::runtime_error(
                 "Molar fraction matrix has wrong dimensions.");
          X result;
          result.context = static_cast<int>(context);
          result.p = p;
          result.T = T;
          auto Sin = S.unchecked<1>();
          auto Cin = C.unchecked<2>();
          for (int iph = 0; iph < X::Model::np; ++iph) {
             result.S[iph] = Sin(iph);
             for (int ic = 0; ic < X::Model::nc; ++ic) {
                result.C[iph][ic] = Cin(iph, ic);
             }
          }
          result.FreeFlow_phase_flowrate.fill(0);
          auto set_outflow = [&]() {
             auto q = py::cast<py::array_t<double, py::array::c_style |
                                                       py::array::forcecast> >(
                 outflow_mass_flowrates);
             if (!(q.ndim() == 1 && q.size() == X::Model::nc))
                throw std::runtime_error(
                    "Mass flowrates vector has wrong dimensions.");
             auto qin = q.unchecked<1>();
             for (int ic = 0; ic < X::Model::nc; ++ic)
                result.FreeFlow_phase_flowrate[ic] = qin(ic);
          };
          switch (context) {
             case Context::gas:
             case Context::liquid:
             case Context::diphasic:
                assert(outflow_mass_flowrates.is_none());
                break;
             case Context::gas_FF_no_liq_outflow:
             case Context::diphasic_FF_no_liq_outflow:
             case Context::diphasic_FF_liq_outflow:
                if (!outflow_mass_flowrates.is_none()) set_outflow();
                break;
             default:
                throw std::runtime_error(
                    "Requested context is not implemented!");
          }
          return result;
       },
       py::arg("context"), py::arg("p"), py::arg("T"), py::arg("S"),
       py::arg("C"), py::arg("outflow_mass_flowrates") = py::none{},
       R"doc(
Construct a state given a specific context and physical parameters.

Parameters
----------

:param context: context (i.e. liquid, gas or diphasic)
:param p: pressure
:param T: temperature
:param S: gaz phase saturations ([Sg, Sl])
:param C: component molar fraction matrix ([[Cga, Cgw], [Cla, Clw]])

)doc");
}
