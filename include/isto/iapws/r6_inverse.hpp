#pragma once
#include "r6.hpp"

/*
    r6.hpp has

      pressure                       (density, temperature)
      massic_internal_energy         (density, temperature)
      massic_entropy                 (density, temperature)
      massic_enthalpy                (density, temperature)
      massic_isochoric_heat_capacity (density, temperature)
      massic_isobaric_heat_capacity  (density, temperature)
      massic_gibbs_free_energy       (density, temperature)
      speed_of_sound                 (density, temperature)

    and

      pressure                       (temperature, density)
      massic_internal_energy         (temperature, density)
      massic_entropy                 (temperature, density)
      massic_enthalpy                (temperature, density)
      massic_isochoric_heat_capacity (temperature, density)
      massic_isobaric_heat_capacity  (temperature, density)
      massic_gibbs_free_energy       (temperature, density)
      speed_of_sound                 (temperature, density)

    Here, we define

      density (pressure, temperature)
      density (temperature, pressure)

      temperature (pressure, density)
      temperature (density, pressure)

    by inverting the firsts.

 */
    namespace
isto::iapws::r6
{
    inline namespace
r6_95_2016
{

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_PRESSURE    pressure_t <T>
#define ISTO_IAPWS_TEMPERATURE temperature_t <T>
#else
#define ISTO_IAPWS_PRESSURE    T
#define ISTO_IAPWS_TEMPERATURE T
    using std::pow;
#endif

    template <class T>
    auto
density (ISTO_IAPWS_PRESSURE const& pressure, ISTO_IAPWS_TEMPERATURE const& temperature)
{
        auto const
    pi = pressure / density / massic_gas_constant / temperature;
        auto const
    tau = critical_temperature / temperature;
        auto const
    delta = zhang ([=](auto delta){ return 1 + delta * detail::phi_r_d (delta, tau) - pi; }, AA, BB, CVG);
    // OR use newton with R7 as a starting point.
    return delta * critical_density;
}
} // inline namespace r6_95_2016
} // namespace isto::iapws::r6
