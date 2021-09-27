#pragma once
#include "r6.hpp"
#include "r7.hpp"
#include <isto/root_finding/root_finding.hpp>
    using namespace isto::root_finding;

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

    Here, we define

      density (pressure, temperature)
      density (temperature, pressure)

    by inverting the above
 */
    namespace
isto::iapws::r6_inverse
{
    inline namespace
r6_95_2016
{
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
density_pt (
      ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_T2 const& temperature
    , ISTO_IAPWS_P3 const& epsilon = 1e-6 ISTO_IAPWS_U_P
){
        using namespace r6;
        auto
    info = info_iterations_t {};
        auto
    tau = critical_temperature / temperature;
        auto&&
    density = newton (
          [=](auto density){ return r6::pressure_dt (density, temperature) - pressure; }
        , [=](auto density)
          {
                using namespace detail;
                auto
            delta = density / critical_density;
            return (1 + 2 * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau)) * massic_gas_constant * temperature;
          }
        , r7::density_pt (pressure, temperature)
        , [=](auto x)
          { 
                using std::abs;
            return abs (x) < epsilon; 
          }
        , {} // options
        , info
    );
    return density;
}
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
density_tp (
      ISTO_IAPWS_T2 const& temperature
    , ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_P3 const& epsilon = 1e-6 ISTO_IAPWS_U_P
){
    return density_pt (pressure, temperature, epsilon);
}
/*
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
temperature_dp (
      ISTO_IAPWS_D2 const& density
    , ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_T3 const& epsilon = 1e-6 ISTO_IAPWS_U_T
){
        using namespace r6;
        auto
    info = info_iterations_t {};
        auto
    delta = density / critical_density;
        auto&&
    temperature = newton (
          [=](auto temperature){ return r6::pressure_dt (density, temperature) - pressure; }
        , [=](auto temperature)
          {
                using namespace detail;
                auto
            tau = critical_temperature / temperature;
            return (1 + delta * phi_r_d (delta, tau) - delta * tau * phi_r_dt (delta, tau)) * massic_gas_constant * density;
          }
        , r7::temperature_dp (density, pressure)
        , [=](auto x)
          { 
            return fabs (x) < epsilon; 
          }
        , {} // options
        , info
    );
    return temperature;
}
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
temperature_pd (
      ISTO_IAPWS_P2 const& pressure
    , ISTO_IAPWS_D1 const& density
    , ISTO_IAPWS_T3 const& epsilon = 1e-6 ISTO_IAPWS_U_T
){
    return temperature_dp (density, pressure, epsilon);
}
*/
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
density (
      ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_T2 const& temperature
    , ISTO_IAPWS_P3 const& epsilon = 1e-6 ISTO_IAPWS_U_P
){
    return density_pt (pressure, temperature, epsilon);
}
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
density (
      ISTO_IAPWS_T2 const& temperature
    , ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_P3 const& epsilon = 1e-6 ISTO_IAPWS_U_P
){
    return density_tp (temperature, pressure, epsilon);
}
/*
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
temperature (
      ISTO_IAPWS_D2 const& density
    , ISTO_IAPWS_P1 const& pressure
    , ISTO_IAPWS_T3 const& epsilon = 1e-6 ISTO_IAPWS_U_T
){
    return temperature_dp (density, pressure, epsilon);
}
    template <class T, class U, class V = std::common_type_t <T, U>>
    auto
temperature (
      ISTO_IAPWS_P2 const& pressure
    , ISTO_IAPWS_D1 const& density
    , ISTO_IAPWS_T3 const& epsilon = 1e-6 ISTO_IAPWS_U_T
){
    return temperature_pd (pressure, density, epsilon);
}
*/
#endif
} // inline namespace r6_95_2016
} // namespace isto::iapws::r6_inverse
