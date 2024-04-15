#pragma once
#include "diagram.hpp"
#include "../include/calculisto/iapws/g12.hpp"
    using namespace calculisto::iapws;
    const auto
topic_g12 = topic_t <double (double, double)>
{
      .name = "g12"
    , .ranges = {
          {
              "ylin"
            , { 181.16, 273.16 }
            , { 0.,     400e6  }
            , { false,  false  }
            , exclusion_tp_supercooled
          }
        , {
              "ylog"
            , { 181.16, 273.16 }
            , { 1.,     400e6  }
            , { false,  true   }
            , exclusion_tp_supercooled
          }
      }
    , .graphs = {
          { "density",                              "temperature", "pressure", g12::density_tp                        <double, double> }
        , { "massic_entropy",                       "temperature", "pressure", g12::massic_entropy_tp                 <double, double>, true }
        , { "isothermal_compressibility",           "temperature", "pressure", g12::isothermal_compressibility_tp     <double, double> }
        , { "thermal_expansion_coefficient",        "temperature", "pressure", g12::thermal_expansion_coefficient_tp  <double, double> }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", g12::massic_isobaric_heat_capacity_tp  <double, double> }
        , { "massic_isochoric_heat_capacity",       "temperature", "pressure", g12::massic_isochoric_heat_capacity_tp <double, double> }
        , { "speed_of_sound",                       "temperature", "pressure", g12::speed_of_sound_tp                 <double, double> }
      }
    , .title = "G12: Guideline on Thermodynamic Properties of Supercooled Water "
    , .description = R"(
<p>Provides direct calculations of the properties of supercooled water in terms of (pressure, temperature).</p>
<p>Range: From the homogeneous ice nucleation temperature to the melting temperature, and from 0 MPa to 400 MPa.</p>
<p>
Implements the following:
<ul>
<li>Vincent Holten, Jan V. Sengers, and Mikhail A. Anisimov. "Equation of State for Supercooled Water at Pressures up to 400 MPa.", J. Phys. Chem. Ref. Data 1 December 2014; 43 (4): 043101. <a href="https://doi.org/10.1063/1.4895593">https://doi.org/10.1063/1.4895593</a></li>
</ul>
</p>
      )"
    , .substance = "Ice Ih"
};
