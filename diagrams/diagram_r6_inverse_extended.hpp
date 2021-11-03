#pragma once
#include "diagram_r6_inverse.hpp"
/*
    auto
exclusion_r6_inverse_extended = exclusion_r6_inverse_t {};
*/
    const auto
topic_r6_inverse_extended = topic_t <double (double, double)>
{
      .name = "r6_inverse_extended"
    , .ranges = {
          {
              "ylin"
            , { 173.16, 5073.16 }
            , { 0.,     100e9   }
            , { false,  false   }
            , exclusion_tp_fluid
          }
        , {
              "ylog"
            , { 173.16, 5073.16 }
            , { 1.    , 100e9   }
            , { false,  true    }
            , exclusion_tp_fluid
          }
      }
    
    , .graphs = {
          { "density",                              "temperature", "pressure", [](double d, double) { return d; }                           }
        , { "massic_enthalpy",                      "temperature", "pressure", r6::massic_enthalpy_dt                      <double, double> }
        , { "massic_internal_energy",               "temperature", "pressure", r6::massic_internal_energy_dt               <double, double> }
        , { "massic_gibbs_free_energy",             "temperature", "pressure", r6::massic_gibbs_free_energy_dt             <double, double> }
        , { "massic_entropy",                       "temperature", "pressure", r6::massic_entropy_dt                       <double, double> }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", r6::massic_isobaric_heat_capacity_dt        <double, double> }
        , { "massic_isochoric_heat_capacity",       "temperature", "pressure", r6::massic_isochoric_heat_capacity_dt       <double, double> }
        , { "speed_of_sound",                       "temperature", "pressure", r6::speed_of_sound_dt                       <double, double> }
        //, { "isobaric_cubic_expansion_coefficient", "temperature", "pressure", r6::isobaric_cubic_expansion_coefficient_dt <double, double> }
        //, { "isothermal_compressibility",           "temperature", "pressure", r6::isothermal_compressibility_dt           <double, double> }
        , { "relative_pressure_coefficient",        "temperature", "pressure", r6::relative_pressure_coefficient_dt        <double, double> }
        , { "isothermal_stress_coefficient",        "temperature", "pressure", r6::isothermal_stress_coefficient_dt        <double, double> }
      }
    , .title = "R6 inverse, extended to 5000 K, 100 GPa"
    , .description = R"(
<p>This is R6 inverse over an extended (P, T) range.</p>
<p>The phase boundary for P &gt; 20 GPa is totally arbitrary, as no model is available in this range (AFAIK).</p>
<p>Range: 173.15 K &le; T &le; 5073.15 K and 0 MPa &le; P &le; 100 GPa </p>
      )"
};
