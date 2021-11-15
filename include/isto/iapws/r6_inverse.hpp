#pragma once
#include "r6.hpp"
#include "r7.hpp"
#include <isto/root_finding/root_finding.hpp>
    using namespace isto::root_finding;

    namespace
isto::iapws::r6_inverse
{
    inline namespace
r6_95_2016
{
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_pt (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_D auto const& initial_guess
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
        using namespace r6;
        auto
    tau = critical_temperature / temperature;
    return newton (
          [=](auto density){ return r6::pressure_dt (density, temperature) - pressure; }
        , [=](auto density)
          {
                using namespace detail;
                auto
            delta = density / critical_density;
            return (1 + 2 * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau)) * massic_gas_constant * temperature;
          }
        , initial_guess
        , [=](auto x)
          { 
                using std::abs;
            return abs (x) < epsilon; 
          }
        , {} // options
        , info
    );
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_pt (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_D auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, 1e-6 ISTO_IAPWS_U_P, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_pt (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, r7::density_pt (pressure, temperature), info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_tp (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_D auto const& initial_guess
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, epsilon, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_tp (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_D auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density_tp (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, info);
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_D auto const& initial_guess
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, epsilon, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_D auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, r7::density_pt (pressure, temperature), epsilon, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_T auto const& temperature
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_D auto const& initial_guess
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
    return density_tp (temperature, pressure, initial_guess, epsilon, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_D auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return density_tp (temperature, pressure, initial_guess, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , ISTO_IAPWS_P auto const& epsilon
    , info_t <InfoTag> info = info::none
){
    return density_tp (temperature, pressure, r7::density_tp (temperature, pressure), epsilon, info);
}
    template <info_tag_t InfoTag = info::tag::none>
    auto
density (
      ISTO_IAPWS_T auto const& temperature
    , ISTO_IAPWS_P auto const& pressure
    , info_t <InfoTag> info = info::none
){
    return density_tp (temperature, pressure, info);
}
#endif
} // inline namespace r6_95_2016
} // namespace isto::iapws::r6_inverse
