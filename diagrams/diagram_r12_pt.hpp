#pragma once
#include "diagram.hpp"
#include "../include/isto/iapws/r12.hpp"
    using namespace isto::iapws;
    const auto
topic_r12_pt = topic_t <double (double, double)>
{
      .name = "r12_pt"
    , .ranges = {
          {
              "ylin"
            , { 173.16, 1173.16 }
            , { 0.,     1000e6  }
            , { false, false  }
            , exclusion_tp_r12
          }
        , {
              "ylog"
            , { 173.16, 1173.16 }
            , { 1.,     1000e6  }
            , { false,  true   }
            , exclusion_tp_r12
          }
      }
    , .graphs = {
          { "viscosity", "temperature", "pressure", r12::viscosity_dt <double, double> }
      }
    , .title = "R12, viscosity in the temperature-pressure plane"
    , .description = R"(
<p>Viscosity in th P-T plane using r6 inverse</p>
<p>Range: 0 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 1000 MPa</p>
<p>
</p>
)"
};
