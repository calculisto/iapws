#pragma once
#include "diagram.hpp"
#include "../include/isto/iapws/r7.hpp"
    using namespace isto::iapws;
    const auto
topic_r7 = topic_t <double (double, double)>
{
      .name = "r7"
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
          { "density",                              "temperature", "pressure", r7::density_tp                              <double, double> }
        , { "massic_enthalpy",                      "temperature", "pressure", r7::massic_enthalpy_tp                      <double, double> }
        , { "massic_internal_energy",               "temperature", "pressure", r7::massic_internal_energy_tp               <double, double> }
        , { "massic_entropy",                       "temperature", "pressure", r7::massic_entropy_tp                       <double, double> }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", r7::massic_isobaric_heat_capacity_tp        <double, double> }
        , { "massic_isochoric_heat_capacity",       "temperature", "pressure", r7::massic_isochoric_heat_capacity_tp       <double, double> }
        , { "speed_of_sound",                       "temperature", "pressure", r7::speed_of_sound_tp                       <double, double> }
        , { "isobaric_cubic_expansion_coefficient", "temperature", "pressure", r7::isobaric_cubic_expansion_coefficient_tp <double, double> }
        , { "isothermal_compressibility",           "temperature", "pressure", r7::isothermal_compressibility_tp           <double, double> }
        , { "relative_pressure_coefficient",        "temperature", "pressure", r7::relative_pressure_coefficient_tp        <double, double> }
        , { "isothermal_stress_coefficient",        "temperature", "pressure", r7::isothermal_stress_coefficient_tp        <double, double> }
      }
    , .title = "R7: The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam"
    , .description = R"(
<p>Provides direct calculations of the properties of ordinary water (i.e. water and stean) in terms of (pressure, temperature).</p>
<p>It is an approximation of R6.</p>
<p>Range: (173.15 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 100 MPa) and (1073.15 K &le; T &le; 2273.15 K and 0 MPa &le; P &le; 50 MPa)</p>
<p>
Implements the following:
<ul>
<li>Wagner, W., Cooper, J. R., Dittmann, A., Kijima, J., Kretzschmar, H., Kruse, A., Mareš, R., Oguchi, K., Sato, H., Stöcker, I., Sǐfner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., and Willkommen, T. (January 1, 2000). "The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam ." ASME. J. Eng. Gas Turbines Power. January 2000; 122(1): 150–184. <a href="https://doi.org/10.1115/1.483186">https://doi.org/10.1115/1.483186</a></li>
<li>Kretzschmar, H., Cooper, J. R., Dittmann, A., Friend, D. G., Gallagher, J. S., Knobloch, K., Mareš, R., Miyagawa, K., Stöcker, I., Trübenbach, J., Wagner, W., and Willkommen, T. (June 22, 2004). "Supplementary Backward Equations for Pressure as a Function of Enthalpy and Entropy li(h,s) to the Industrial Formulation IAPWS-IF97 for Water and Steam." ASME. J. Eng. Gas Turbines Power. July 2006; 128(3): 702–713. <a href="https://doi.org/10.1115/1.1915392">https://doi.org/10.1115/1.1915392</a></li>
<li>Kretzschmar, H., Cooper, J. R., Dittmann, A., Friend, D. G., Gallagher, J. S., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Stöcker, I., Wagner, W., and Weber, I. (January 10, 2006). "Supplementary Backward Equations T(li,h)⁠, v(li,h)⁠, and T(li,s)⁠, v(li,s) for the Critical and Supercritical Regions (Region 3) of the Industrial Formulation IAPWS-IF97 for Water and Steam." ASME. J. Eng. Gas Turbines Power. January 2007; 129(1): 294–303. <a href="https://doi.org/10.1115/1.2181598">https://doi.org/10.1115/1.2181598</a></li>
<li>Kretzschmar, H., Cooper, J. R., Gallagher, J. S., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Span, R., Stöcker, I., Wagner, W., and Weber, I. (January 16, 2007). "Supplementary Backward Equations li(h,s) for the Critical and Supercritical Regions (Region 3), and Equations for the Two-Phase Region and Region Boundaries of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam." ASME. J. Eng. Gas Turbines Power. October 2007; 129(4): 1125–1137. <a href="https://doi.org/10.1115/1.2719267">https://doi.org/10.1115/1.2719267</a></li>
<li>Kretzschmar, H., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Span, R., Stöcker, I., Wagner, W., and Weber, I. (April 13, 2009). "Supplementary Backward Equations v(li,T) for the Critical and Supercritical Regions (Region 3) of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam." ASME. J. Eng. Gas Turbines Power. July 2009; 131(4): 043101. <a href="https://doi.org/10.1115/1.3028630">https://doi.org/10.1115/1.3028630</a></li>
</ul>
</p>
      )"
};
