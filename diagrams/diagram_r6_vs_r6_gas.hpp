#include "diagram.hpp"
#include <optional>
#include "../include/calculisto/iapws/r6.hpp"
    struct
delta2
{
        std::function <double (double, double)>
    f6;
        std::function <double (double, double)>
    f6_gas;
        template <class F6, class F6_gas>
    delta2 (F6 a, F6_gas b)
        : f6 (a)
        , f6_gas (b)
    {}
        auto
    operator () (double d, double t)
    {
            const auto
        a = f6 (d, t);
            const auto
        b = f6_gas (d, t);
        return fabs ((a - b) / a);
    }
};

    const auto
topic_r6_vs_r6_gas = topic_t <double (double, double)>
{
      .name = "r6_vs_r6_gas"
    , .ranges = {
          {
              "ylin"
            , { 0.,     55.   }
            , { 273.16, 2273.16 }
            , { false,  false   }
            , exclusion_dt
          }
      }
    
    , .graphs = {
          { "pressure",                             "density", "temperature", delta2 { r6::pressure_dt                             <double, double>, r6_gas::pressure_dt                             <double, double> }}
        , { "massic_enthalpy",                      "density", "temperature", delta2 { r6::massic_enthalpy_dt                      <double, double>, r6_gas::massic_enthalpy_dt                      <double, double> }}
        , { "massic_internal_energy",               "density", "temperature", delta2 { r6::massic_internal_energy_dt               <double, double>, r6_gas::massic_internal_energy_dt               <double, double> }}
        , { "massic_gibbs_free_energy",             "density", "temperature", delta2 { r6::massic_gibbs_free_energy_dt             <double, double>, r6_gas::massic_gibbs_free_energy_dt             <double, double> }}
        , { "massic_entropy",                       "density", "temperature", delta2 { r6::massic_entropy_dt                       <double, double>, r6_gas::massic_entropy_dt                       <double, double> }}
        , { "massic_isobaric_heat_capacity",        "density", "temperature", delta2 { r6::massic_isobaric_heat_capacity_dt        <double, double>, r6_gas::massic_isobaric_heat_capacity_dt        <double, double> }}
        , { "massic_isochoric_heat_capacity",       "density", "temperature", delta2 { r6::massic_isochoric_heat_capacity_dt       <double, double>, r6_gas::massic_isochoric_heat_capacity_dt       <double, double> }}
        , { "speed_of_sound",                       "density", "temperature", delta2 { r6::speed_of_sound_dt                       <double, double>, r6_gas::speed_of_sound_dt                       <double, double> }}
        , { "relative_pressure_coefficient",        "density", "temperature", delta2 { r6::relative_pressure_coefficient_dt        <double, double>, r6_gas::relative_pressure_coefficient_dt        <double, double> }}
        , { "isothermal_stress_coefficient",        "density", "temperature", delta2 { r6::isothermal_stress_coefficient_dt        <double, double>, r6_gas::isothermal_stress_coefficient_dt        <double, double> }}
        //, { "isobaric_cubic_expansion_coefficient", "density", "temperature", r6::isobaric_cubic_expansion_coefficient_dt <double, double> }
        //, { "isothermal_compressibility",           "density", "temperature", r6::isothermal_compressibility_dt           <double, double> }
      }
    , .title = "R6 vs. R6_gas: Compare R6 gas equation to the generic fluid equation"
    , .description = R"(
<p></p>
<p>Range : 273 K &le; T &lt; 2273 K and 0 kg/m<sup>3</sup> &le; &rho; &le; 55 kg/m<sup>3</sup></p>
<p>
Implementsn:
<ul>
<li>W. Wagner and A. Pruß , "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use", Journal of Physical and Chemical Reference Data 31, 387-535 (2002) <a href="https://doi.org/10.1063/1.1461829">https://doi.org/10.1063/1.1461829</a></li>
</ul>
</p>
      )"
    , .substance = "Water"
};

