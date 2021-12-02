#include "diagram.hpp"
#include <optional>
#include "../include/isto/iapws/r6.hpp"
    const auto
topic_r6 = topic_t <double (double, double)>
{
      .name = "r6"
    , .ranges = {
          {
              "ylin"
            , { 0.,     1000.   }
            , { 273.16, 2273.16 }
            , { false,  false   }
            , exclusion_dt
          }
      }
    
    , .graphs = {
          { "pressure",                             "density", "temperature", r6::pressure_dt                             <double, double> }
        , { "massic_enthalpy",                      "density", "temperature", r6::massic_enthalpy_dt                      <double, double> }
        , { "massic_internal_energy",               "density", "temperature", r6::massic_internal_energy_dt               <double, double> }
        , { "massic_gibbs_free_energy",             "density", "temperature", r6::massic_gibbs_free_energy_dt             <double, double> }
        , { "massic_entropy",                       "density", "temperature", r6::massic_entropy_dt                       <double, double> }
        , { "massic_isobaric_heat_capacity",        "density", "temperature", r6::massic_isobaric_heat_capacity_dt        <double, double> }
        , { "massic_isochoric_heat_capacity",       "density", "temperature", r6::massic_isochoric_heat_capacity_dt       <double, double> }
        , { "speed_of_sound",                       "density", "temperature", r6::speed_of_sound_dt                       <double, double> }
        //, { "isobaric_cubic_expansion_coefficient", "density", "temperature", r6::isobaric_cubic_expansion_coefficient_dt <double, double> }
        //, { "isothermal_compressibility",           "density", "temperature", r6::isothermal_compressibility_dt           <double, double> }
        , { "relative_pressure_coefficient",        "density", "temperature", r6::relative_pressure_coefficient_dt        <double, double> }
        , { "isothermal_stress_coefficient",        "density", "temperature", r6::isothermal_stress_coefficient_dt        <double, double> }
      }
    , .title = "R6: The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use"
    , .description = R"(
<p>Provides direct calculations of the properties of ordinary water (i.e. water and stean) in terms of (density, temperature).</p>
<p>Range : 273.16 K &le; T &lt; 2273.16 K and 0 kg/m<sup>3</sup> &le; &rho; &le; 1000 kg/m<sup>3</sup></p>
<p>
Implements the following:
<ul>
<li>W. Wagner and A. Pruß , "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use", Journal of Physical and Chemical Reference Data 31, 387-535 (2002) <a href="https://doi.org/10.1063/1.1461829">https://doi.org/10.1063/1.1461829</a></li>
</ul>
</p>
      )"
    , .substance = "Water"
};

