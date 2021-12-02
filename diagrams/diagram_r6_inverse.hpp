#pragma once
#include "diagram.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws;
    const auto
topic_r6_inverse = topic_t <double (double, double)>
{
      .name = "r6_inverse"
    , .ranges = {
          {
              "ylin"
            , { 173.16, 1273.16 }
            , { 0.,     1000e6  }
            , { false,  false   }
            , exclusion_tp_fluid
          }
        , {
              "ylog"
            , { 173.16, 2273.16 }
            , { 1.    , 1000e6  }
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
    , .title = "R6 inverse"
    , .description = R"(
<p>Provides indirect calculations of the properties of ordinary water in terms 
of (pressure, temperature) by inverting the P(&rho;, T) relation of R6.</p>
<p>Range: 173.15 K &le; T &le; 1273.15 K and 0 MPa &le; P &le; 1000 MPa </p>
    )"
    , .substance = "Water"
};

