#pragma once
#include "detail/common.hpp"
#include <complex>
#include "detail/template_pow.hpp"

    namespace 
isto::iapws::r10
{
    inline namespace 
r10_06_2009
{
    namespace
detail
{
    constexpr auto
triple_point_pressure = 611.657 ISTO_IAPWS_U_P;
    constexpr auto
normal_pressure = 101325 ISTO_IAPWS_U_P;
    constexpr auto
triple_point_temperature = 273.16 ISTO_IAPWS_U_T;

    constexpr auto
g00 = -0.632020233335886e6   ISTO_IAPWS_U_H;
   constexpr auto
g01 =  0.655022213658955     ISTO_IAPWS_U_H;
   constexpr auto
g02 = -0.189369929326131e-7  ISTO_IAPWS_U_H;
   constexpr auto
g03 =  0.339746123271053e-14 ISTO_IAPWS_U_H;
   constexpr auto
g04 = -0.556464869058991e-21 ISTO_IAPWS_U_H;

    constexpr auto
s0_abs = 0.18913e3 ISTO_IAPWS_U_S;
    constexpr auto
s0 = -0.332733756492168e4 ISTO_IAPWS_U_S;

    constexpr auto
t1  = std::complex {  0.368017112855051e-1,   0.510878114959572e-1  };
   constexpr auto
r1  = std::complex {  0.447050716285388e2,    0.656876847463481e2   } ISTO_IAPWS_U_S;
   constexpr auto
t2  = std::complex {  0.337315741065416,      0.335449415919309     };
    constexpr auto
r20 = std::complex { -0.725974574329220e2,   -0.781008427112870e2   } ISTO_IAPWS_U_S;
    constexpr auto
r21 = std::complex { -0.557107698030123e-4,   0.464578634580806e-4  } ISTO_IAPWS_U_S;
    constexpr auto
r22 = std::complex {  0.234801409215913e-10, -0.285651142904972e-10 } ISTO_IAPWS_U_S;

    constexpr auto
pi0 = 101325. / 611.657;

    constexpr auto
g0 (auto const& pi)
{
        using std::pow;
    return
          g00
        + g01 * pow (pi - pi0, 1)
        + g02 * pow (pi - pi0, 2)
        + g03 * pow (pi - pi0, 3)
        + g04 * pow (pi - pi0, 4)
    ;
}

    constexpr auto
g0_p (auto const& pi)
{
        using std::pow;
    return
          g01 * 1 / triple_point_pressure
        + g02 * 2 / triple_point_pressure * pow (pi - pi0, 1)
        + g03 * 3 / triple_point_pressure * pow (pi - pi0, 2)
        + g04 * 4 / triple_point_pressure * pow (pi - pi0, 3)
    ;
}

    constexpr auto
g0_pp (auto const& pi)
{
        using std::pow;
    return
        + g02 * 2 * 1 / triple_point_pressure / triple_point_pressure
        + g03 * 3 * 2 / triple_point_pressure / triple_point_pressure * pow (pi - pi0, 1)
        + g04 * 4 * 3 / triple_point_pressure / triple_point_pressure * pow (pi - pi0, 2)
    ;
}

    constexpr auto
r2 (auto const& pi)
{
        using std::pow;
    return
          r20
        + r21 * pow (pi - pi0, 1)
        + r22 * pow (pi - pi0, 2)
    ;
}

    constexpr auto
r2_p (auto const& pi)
{
        using std::pow;
    return
          r21 * 1. / triple_point_pressure * pow (pi - pi0, 0)
        + r22 * 2. / triple_point_pressure * pow (pi - pi0, 1)
    ;
}

    constexpr auto
r2_pp = r22 * 2. / triple_point_pressure / triple_point_pressure;

    constexpr auto
g (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return g0 (pi) - s0 * temperature + triple_point_temperature * real (
          r1      * ((t1 - tau) * log (t1 - tau) + (t1 + tau) * log (t1 + tau) - 2. * t1 * log (t1) - tau * tau / t1)
        + r2 (pi) * ((t2 - tau) * log (t2 - tau) + (t2 + tau) * log (t2 + tau) - 2. * t2 * log (t2) - tau * tau / t2)
    );
}

    constexpr auto
g_t (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return - s0 + real (
          r1      * (-log (t1 - tau) + log (t1 + tau) - 2. * tau / t1)
        + r2 (pi) * (-log (t2 - tau) + log (t2 + tau) - 2. * tau / t2)
    );
}

    constexpr auto
g_p (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return g0_p (pi) + triple_point_temperature * real (
        + r2_p (pi) * ((t2 - tau) * log (t2 - tau) + (t2 + tau) * log (t2 + tau) - 2. * t2 * log (t2) - tau * tau / t2)
    );
}

    constexpr auto
g_tt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return 1. / triple_point_temperature * real (
          r1      * (1. / (t1 - tau) + 1. / (t1 + tau) - 2. / t1)
        + r2 (pi) * (1. / (t2 - tau) + 1. / (t2 + tau) - 2. / t2)
    );
}

    constexpr auto
g_tp (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return real (
        + r2_p (pi) * (-log (t2 - tau) + log (t2 + tau) - 2. * tau / t2)
    );
}

    constexpr auto
g_pp (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using std::real, std::log;
        auto const
    tau = temperature / triple_point_temperature;
        auto const
    pi = pressure / triple_point_pressure;
    return g0_pp (pi) + triple_point_temperature * real (
        + r2_pp * ((t2 - tau) * log (t2 - tau) + (t2 + tau) * log (t2 + tau) - 2. * t2 * log (t2) - tau * tau / t2)
    );
}
} // namespace detail
    using namespace detail;
    constexpr auto
massic_volume_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return g_p (pressure, temperature);
}
    constexpr auto
massic_volume_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_volume_pt (pressure, temperature);
}
    constexpr auto
density_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return 1 / massic_volume_pt (pressure, temperature);
}
    constexpr auto
density_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return density_pt (pressure, temperature);
}
    constexpr auto
massic_entropy_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return -g_t (pressure, temperature);
}
    constexpr auto
massic_entropy_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_entropy_pt (pressure, temperature);
}
    constexpr auto
massic_isobaric_heat_capacity_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return -g_tt (pressure, temperature) * temperature;
}
    constexpr auto
massic_isobaric_heat_capacity_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    constexpr auto
massic_enthalpy_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return g (pressure, temperature) - temperature * g_t (pressure, temperature);
}
    constexpr auto
massic_enthalpy_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    constexpr auto
massic_internal_energy_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return g (pressure, temperature) 
        - temperature * g_t (pressure, temperature)
        - pressure * g_p (pressure, temperature)
    ;
}
    constexpr auto
massic_internal_energy_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    constexpr auto
massic_helmholtz_energy_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return g (pressure, temperature) 
        - pressure * g_p (pressure, temperature)
    ;
}
    constexpr auto
massic_helmholtz_energy_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_helmholtz_energy_pt (pressure, temperature);
}
    constexpr auto
cubic_expansion_coefficient_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return g_tp (pressure, temperature) / g_p (pressure, temperature);
}
    constexpr auto
cubic_expansion_coefficient_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return cubic_expansion_coefficient_pt (pressure, temperature);
}
    constexpr auto
pressure_coefficient_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return -g_tp (pressure, temperature) / g_pp (pressure, temperature);
}
    constexpr auto
pressure_coefficient_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return pressure_coefficient_pt (pressure, temperature);
}
    constexpr auto
isothermal_compressibility_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return -g_pp (pressure, temperature) / g_p (pressure, temperature);
}
    constexpr auto
isothermal_compressibility_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isentropic_compressibility_pt (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
        using isto::template_pow::pow;
    return (pow <2> (g_tp (pressure, temperature)) - g_tt (pressure, temperature) * g_pp (pressure, temperature)) / (g_p (pressure, temperature) * g_tt (pressure, temperature));
}
    constexpr auto
isentropic_compressibility_tp (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return isentropic_compressibility_pt (pressure, temperature);
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    using namespace detail;
    constexpr auto
massic_volume (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_volume_pt (pressure, temperature);
}
    constexpr auto
massic_volume (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_volume_pt (pressure, temperature);
}
    constexpr auto
density (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return density_pt (pressure, temperature);
}
    constexpr auto
density (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return density_pt (pressure, temperature);
}
    constexpr auto
massic_entropy (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_entropy_pt (pressure, temperature);
}
    constexpr auto
massic_entropy (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_entropy_pt (pressure, temperature);
}
    constexpr auto
massic_isobaric_heat_capacity (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    constexpr auto
massic_isobaric_heat_capacity (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    constexpr auto
massic_enthalpy (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    constexpr auto
massic_enthalpy (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    constexpr auto
massic_internal_energy (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    constexpr auto
massic_internal_energy (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    constexpr auto
massic_helmholtz_energy (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return massic_helmholtz_energy_pt (pressure, temperature);
}
    constexpr auto
massic_helmholtz_energy (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return massic_helmholtz_energy_pt (pressure, temperature);
}
    constexpr auto
cubic_expansion_coefficient (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return cubic_expansion_coefficient_pt (pressure, temperature);
}
    constexpr auto
cubic_expansion_coefficient (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return cubic_expansion_coefficient_pt (pressure, temperature);
}
    constexpr auto
pressure_coefficient (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return pressure_coefficient_pt (pressure, temperature);
}
    constexpr auto
pressure_coefficient (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return pressure_coefficient_pt (pressure, temperature);
}
    constexpr auto
isothermal_compressibility (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isothermal_compressibility (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isentropic_compressibility (ISTO_IAPWS_P auto const& pressure, ISTO_IAPWS_T auto const& temperature)
{
    return isentropic_compressibility_pt (pressure, temperature);
}
    constexpr auto
isentropic_compressibility (ISTO_IAPWS_T auto const& temperature, ISTO_IAPWS_P auto const& pressure)
{
    return isentropic_compressibility_pt (pressure, temperature);
}
#endif
} // namespace r10_06_2009
} // namespace isto::iapws::r10
