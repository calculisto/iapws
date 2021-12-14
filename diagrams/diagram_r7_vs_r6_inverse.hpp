#include "diagram.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws;

    double
ident1 (double d, double)
{
    return d;
}

    struct
delta
{
        std::function <double (double, double)>
    f6;
        std::function <double (double, double)>
    f7;
        template <class F6, class F7>
    delta (F6 a, F7 b)
        : f6 (a)
        , f7 (b)
    {}
        auto
    operator () (double t, double p, double d)
    {
            const auto
        a = f6 (d, t);
            const auto
        b = f7 (t, p);
        return fabs ((a - b) / a);
    }
};


    const auto
topic_r7_vs_r6_inverse = topic_t <double (double, double, double)>
{
      .name = "r7_vs_r6_inverse"
    , .ranges = {
          {
              "ylin"
            , { 273.16, 2273.16 }
            , { 0.,     100e6   }
            , { false,  false   }
            , exclusion_tp_r7
          }
        , {
              "ylog"
            , { 273.16, 2273.16 }
            , { 1.    , 100e6   }
            , { false, true    }
            , exclusion_tp_r7
          }
      }
    
    , .graphs = {
          { "delta_density",                              "temperature", "pressure", delta { ident1,                                                        r7::density_tp                              <double, double> }}
        , { "delta_massic_enthalpy",                      "temperature", "pressure", delta { r6::massic_enthalpy_dt                      <double, double>, r7::massic_enthalpy_tp                      <double, double> }}
        , { "delta_massic_internal_energy",               "temperature", "pressure", delta { r6::massic_internal_energy_dt               <double, double>, r7::massic_internal_energy_tp               <double, double> }}
        , { "delta_massic_entropy",                       "temperature", "pressure", delta { r6::massic_entropy_dt                       <double, double>, r7::massic_entropy_tp                       <double, double> }}
        , { "delta_massic_isobaric_heat_capacity",        "temperature", "pressure", delta { r6::massic_isobaric_heat_capacity_dt        <double, double>, r7::massic_isobaric_heat_capacity_tp        <double, double> }}
        , { "delta_massic_isochoric_heat_capacity",       "temperature", "pressure", delta { r6::massic_isochoric_heat_capacity_dt       <double, double>, r7::massic_isochoric_heat_capacity_tp       <double, double> }}
        , { "delta_speed_of_sound",                       "temperature", "pressure", delta { r6::speed_of_sound_dt                       <double, double>, r7::speed_of_sound_tp                       <double, double> }}
        //, { "delta_isobaric_cubic_expansion_coefficient", "temperature", "pressure", delta { r6::isobaric_cubic_expansion_coefficient_dt <double, double>, r7::isobaric_cubic_expansion_coefficient_tp <double, double> }}
        //, { "delta_isothermal_compressibility",           "temperature", "pressure", delta { r6::isothermal_compressibility_dt           <double, double>, r7::isothermal_compressibility_tp           <double, double> }}
        , { "delta_relative_pressure_coefficient",        "temperature", "pressure", delta { r6::relative_pressure_coefficient_dt        <double, double>, r7::relative_pressure_coefficient_tp        <double, double> }}
        , { "delta_isothermal_stress_coefficient",        "temperature", "pressure", delta { r6::isothermal_stress_coefficient_dt        <double, double>, r7::isothermal_stress_coefficient_tp        <double, double> }}
      }
    , .title = "R7 vs. R6 inverse"
    , .description = R"(
<p>Comparison between the values obtained via R7 and R6 inverse.</p>
<p>Range: (173.15 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 100 MPa) and (1073.15 K &le; T &le; 2273.15 K and 0 MPa &le; P &le; 50 MPa)</p>
    )"
    , .substance = "Water"
};


