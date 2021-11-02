#pragma once
#include "diagram.hpp"
#include "../include/isto/iapws/r10.hpp"
#include "../include/isto/iapws/r14.hpp"
    using namespace isto::iapws;

    struct
exclusion_r10_t
    : exclusion_base_t
{
        double
    pmih;
        double
    ps;
        void
    init (double t) override
    {
        pmih  = r14::ih::melting_pressure_t  (t);
        ps    = r14::sublimation_pressure_t  (t);
    };
        bool
    operator () (double t, double p) const override
    {
        return !(
               (p <  611.657   && t <= 273.16 && p > ps)
            || (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)
        );
    };
};
    auto
exclusion_r10 = exclusion_r10_t {};

    const auto
topic_r10 = topic_t <double (double, double)>
{
      .name = "r10"
    , .ranges = {
          {
              "ylin"
            , { 0.,    273.16 }
            , { 0.,    210e6  }
            , { false, false  }
            , exclusion_r10
          }
        , {
              "ylog"
            , { 0.,     273.16 }
            , { 1.,     210e6  }
            , { false,  true   }
            , exclusion_r10
          }
      }
    , .graphs = {
          { "density",                              "temperature", "pressure", r10::density_tp                              <double, double> }
        , { "massic_enthalpy",                      "temperature", "pressure", r10::massic_enthalpy_tp                      <double, double>, true }
        , { "massic_internal_energy",               "temperature", "pressure", r10::massic_internal_energy_tp               <double, double>, true }
        , { "massic_entropy",                       "temperature", "pressure", r10::massic_entropy_tp                       <double, double>, true }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", r10::massic_isobaric_heat_capacity_tp        <double, double> }
        //, { "massic_isochoric_heat_capacity",       "temperature", "pressure", r10::massic_isochoric_heat_capacity_tp       <double, double> }
        //, { "speed_of_sound",                       "temperature", "pressure", r10::speed_of_sound_tp                       <double, double> }
        //, { "isobaric_cubic_expansion_coefficient", "temperature", "pressure", r10::isobaric_cubic_expansion_coefficient_tp <double, double> }
        , { "isothermal_compressibility",           "temperature", "pressure", r10::isothermal_compressibility_tp           <double, double> }
        //, { "relative_pressure_coefficient",        "temperature", "pressure", r10::relative_pressure_coefficient_tp        <double, double> }
        //, { "isothermal_stress_coefficient",        "temperature", "pressure", r10::isothermal_stress_coefficient_tp        <double, double> }
      }
    , .title = "R10: Revised Release on the Equation of State 2006 for H<sub>2</sub>O Ice Ih"
    , .description = R"(
<p>Provides direct calculations of the properties of ordinary water in its solid hexagonal phase I in terms of (pressure, temperature).</p>
<p>Range: 0 K &le; T &lt; 273.16 K and 0 MPa &le; P &le; 210 MPa</p>
<p>
Implements the following:
<ul>
<li>Rainer Feistel, and Wolfgang Wagner. "A New Equation of State for H2O Ice Ih", Journal of Physical and Chemical Reference Data 35, 1021-1047 (2006) <a href="https://doi.org/10.1063/1.2183324">https://doi.org/10.1063/1.2183324</a></li>
</ul>
</p>
      )"
};
