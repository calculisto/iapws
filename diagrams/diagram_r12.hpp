#include "diagram.hpp"
#include "../include/isto/iapws/r12.hpp"
    using namespace isto::iapws;
    const auto
topic_r12 = topic_t <double (double, double)>
{
      .name = "r12"
    , .ranges = {
          {
              "ylin"
            , { 0.,     1000.   }
            , { 273.16, 2173.16 }
            , { false,  false   }
            , exclusion_dt
          }
      }
    
    , .graphs = {
          { "viscosity", "density", "temperature", r12::viscosity_dt <double, double> }
      }
    , .title = "R12: Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance"
    , .description = R"(
<p>Viscosity of liquid and gaseous water in the density-temperature plane.</p>
<p>Range : 273.16 K &le; T &lt; 2273.16 K and 0 kg/m<sup>3</sup> &le; &rho; &le; 1000 kg/m<sup>3</sup></p>
<p>
Implements the following:
<ul>
<li>M. L. Huber, R. A. Perkins, A. Laesecke, D. G. Friend, J. V. Sengers, M. J. Assael I. N. Metaxa, E. Vogel, R. Mareš, and K. Miyagawa , "New International Formulation for the Viscosity of H2O", Journal of Physical and Chemical Reference Data 38, 101-125 (2009) <a href="https://doi.org/10.1063/1.3088050">https://doi.org/10.1063/1.3088050</a></li>
</ul>
</p>
      )"
    , .substance = "Water"
};

