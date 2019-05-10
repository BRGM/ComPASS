#include "Model_wrappers.h"

#include <array>

#include <pybind11/numpy.h>

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;
static_assert(NP == 2, "Wrong numpber of phases.");
static_assert(NC == 1, "Wrong numpber of components.");
static_assert(ComPASS_NUMBER_OF_CONTEXTS==3, "Wrong number of contexts.");
// FIXME: assuming liquid phase is the latest phase
constexpr int GAS_PHASE = 0;
constexpr int LIQUID_PHASE = 1;

enum struct Phase {
    gas = ComPASS_GAS_PHASE,
    liquid = ComPASS_LIQUID_PHASE
};

enum struct Context {
    gas = ComPASS_GAS_CONTEXT,
    liquid = ComPASS_LIQUID_CONTEXT,
    diphasic = ComPASS_DIPHASIC_CONTEXT
};

// Fortran functions
extern "C"
{
    void FluidThermodynamics_molar_density(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_molar_enthalpy(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_dynamic_viscosity(int, double, double, const double *, const double *, double&, double&, double&, double *, double *);
    void FluidThermodynamics_Psat(double, double&, double&);
    void FluidThermodynamics_Tsat(double, double&, double&);
}


void init_model() {}

void finalize_model() {}

template <int PHASE>
inline auto saturate_phase()
{
    auto S = std::array<double, NP>{}; // zero initialization
    S[PHASE] = 1;
    assert(std::accumulate(begin(S), end(S), 0) == 1);
    return S;
}

template <int PHASE>
inline double phase_molar_density(double p, double T)
{
    double xsi, dxsidp, dxsidT;
    double C[NC] = { 1 };
    auto S = saturate_phase<PHASE>();
    double dxsidC[NC] = { 0 };
    auto dxsidS = std::array<double, NP>{}; // zero initialization
    FluidThermodynamics_molar_density(PHASE+1, p, T, C, S.data(), xsi, dxsidp, dxsidT, dxsidC, dxsidS.data());
    return xsi;
}

inline double gas_molar_density(double p, double T)
{
    return phase_molar_density<GAS_PHASE>(p, T);
}

inline double liquid_molar_density(double p, double T)
{
    return phase_molar_density<LIQUID_PHASE>(p, T);
}

template <int PHASE>
inline double phase_molar_enthalpy(double p, double T)
{
    double h, dhdp, dhdT;
    double C[NC] = { 1 };
    auto S = saturate_phase<PHASE>();
    double dhdC[NC] = { 0 };
    auto dhdS = std::array<double, NP>{}; // zero initialization
    FluidThermodynamics_molar_enthalpy(PHASE+1, p, T, C, S.data(), h, dhdp, dhdT, dhdC, dhdS.data());
    return h;
}

inline double gas_molar_enthalpy(double p, double T)
{
    return phase_molar_enthalpy<GAS_PHASE>(p, T);
}

inline double liquid_molar_enthalpy(double p, double T)
{
    return phase_molar_enthalpy<LIQUID_PHASE>(p, T);
}

template <int PHASE>
inline double phase_dynamic_viscosity(double p, double T)
{
    double mu, dmudp, dmudT;
    double C[NC] = { 1 };
    auto S = saturate_phase<PHASE>();
    double dmudC[NC] = { 0 };
    auto dmudS = std::array<double, NP>{}; // zero initialization
    FluidThermodynamics_dynamic_viscosity(PHASE+1, p, T, C, S.data(), mu, dmudp, dmudT, dmudC, dmudS.data());
    return mu;
}

inline double gas_dynamic_viscosity(double p, double T)
{
    return phase_dynamic_viscosity<GAS_PHASE>(p, T);
}

inline double liquid_dynamic_viscosity(double p, double T)
{
    return phase_dynamic_viscosity<LIQUID_PHASE>(p, T);
}

inline double Psat(double T)
{
    double result;
    double dPsatdT;
    FluidThermodynamics_Psat(T, result, dPsatdT);
    return result;
}


inline double Tsat(double p)
{
    double result;
    double dTsatdp;
    FluidThermodynamics_Tsat(p, result, dTsatdp);
    return result;
}

void add_specific_model_wrappers(py::module& module)
{

    module.def("Psat", py::vectorize(Psat));
    module.def("Tsat", py::vectorize(Tsat));
    module.def("gas_molar_density", py::vectorize(gas_molar_density));
    module.def("gas_molar_enthalpy", py::vectorize(gas_molar_enthalpy));
    module.def("gas_dynamic_viscosity", py::vectorize(gas_dynamic_viscosity));
    module.def("liquid_molar_density", py::vectorize(liquid_molar_density));
    module.def("liquid_molar_enthalpy", py::vectorize(liquid_molar_enthalpy));
    module.def("liquid_dynamic_viscosity", py::vectorize(liquid_dynamic_viscosity));
    
    py::enum_<Context>(module, "Context")
        .value("gas", Context::gas)
        .value("liquid", Context::liquid)
        .value("diphasic", Context::diphasic)
    ;

    py::enum_<Phase>(module, "Phase")
        .value("gas", Phase::gas)
        .value("liquid", Phase::liquid)
    ;

    use_context_locker(module);

}
