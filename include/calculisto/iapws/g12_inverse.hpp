#pragma once
#include "g12.hpp"
#include <calculisto/root_finding/root_finding.hpp>
    using namespace calculisto::root_finding;

    namespace
calculisto::iapws::g12_inverse
{
    namespace
detail
{
    constexpr auto
dvdp_tp (auto const& temperature, auto const& pressure)
{
        using namespace g12::detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
        const auto
    omega = 2 + omega_0 * pi;
        const auto
    chi = 1. / (2. / (1 - pow (phi, 2)) - omega);
    return -(
          (tau + 1.) / 2. 
          * (chi * pow (L_p (tau, pi) - omega_0 * phi, 2.) - L_pp (tau, pi) * (phi + 1.))
        - Psi_r_pp (tau, pi)
    ) / rho_0 / rho_0 / R / T_LL;
}

    template <class T>
    constexpr auto
initial_pressure_vt (
      [[maybe_unused]] T const& massic_volume
    , [[maybe_unused]] T const& temperature
){
        const auto
    d = 1. / massic_volume;
    // Use a linear fit for the initial guess.
        const auto
    r = (d - 0.18884 * temperature - 944.073) / 5.01387e-07;
    return r;
}
} // namespace detail
    inline namespace
g12_15
{
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
pressure_vt (
      auto const& massic_volume
    , auto const& temperature
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return newton (
          [=](auto pressure)
          { 
                return g12::massic_volume_tp (temperature, pressure) - massic_volume; 
          }
        , [=](auto pressure)
          {
            return detail::dvdp_tp (temperature, pressure) ;
          }
        , initial_guess
        , { .converged = [](auto curr, auto prev, auto f)
            { 
                return fabs (f) < 1e-8 || fabs ((curr - prev) / curr) < 1e-8; 
            } 
          } // options
        , info
    );
}
    template <info_tag_t InfoTag>
    constexpr auto
pressure_vt (
      auto const& massic_volume
    , auto const& temperature
    , info_t <InfoTag> info = info::none
){
    return pressure_vt (
          massic_volume
        , temperature
        , detail::initial_pressure_vt (massic_volume, temperature)
        , info
    );
}
    constexpr auto
pressure_vt (
      auto const& massic_volume
    , auto const& temperature
){
    return pressure_vt (
          massic_volume
        , temperature
        , detail::initial_pressure_vt (massic_volume, temperature)
        , info::none
    );
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
pressure_dt (
      auto const& density
    , auto const& temperature
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return pressure_vt (1. / density, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag>
    constexpr auto
pressure_dt (
      auto const& density
    , auto const& temperature
    , info_t <InfoTag> info = info::none
){
    return pressure_vt (1. / density, temperature, info);
}
    constexpr auto
pressure_dt (
      auto const& density
    , auto const& temperature
){
    return pressure_vt (1. / density, temperature);
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
pressure_tv (
      auto const& temperature
    , auto const& massic_volume
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return pressure_vt (massic_volume, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag>
    constexpr auto
pressure_tv (
      auto const& temperature
    , auto const& massic_volume
    , info_t <InfoTag> info = info::none
){
    return pressure_vt (massic_volume, temperature, info);
}
    constexpr auto
pressure_tv (
      auto const& temperature
    , auto const& massic_volume
){
    return pressure_vt (massic_volume, temperature);
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
pressure_td (
      auto const& temperature
    , auto const& density
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return pressure_dt (density, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag>
    constexpr auto
pressure_td (
      auto const& temperature
    , auto const& density
    , info_t <InfoTag> info = info::none
){
    return pressure_dt (density, temperature, info);
}
    constexpr auto
pressure_td (
      auto const& temperature
    , auto const& density
){
    return pressure_dt (density, temperature);
}


#define ISTO_IAPWS_G12_INVERSE_GEN(NAME)                                         \
    template <info_tag_t InfoTag = info::tag::none>                              \
    constexpr auto                                                               \
NAME##_vt (                                                                      \
      auto const& massic_volume                                                  \
    , auto const& temperature                                                    \
    , auto const& pressure_initial_guess                                         \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
        const auto                                                               \
    p = pressure_vt (massic_volume, temperature, pressure_initial_guess, info);; \
    if constexpr (InfoTag != info::tag::none)                                    \
    {                                                                            \
        return std::pair { g12::NAME##_pt (p.first, temperature), p.second };    \
    }                                                                            \
    else                                                                         \
    {                                                                            \
        return g12::NAME##_pt (p, temperature);                                  \
    }                                                                            \
}                                                                                \
    template <info_tag_t InfoTag>                                                \
    constexpr auto                                                               \
NAME##_vt (                                                                      \
      auto const& massic_volume                                                  \
    , auto const& temperature                                                    \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_vt (                                                           \
          massic_volume                                                          \
        , temperature                                                            \
        , detail::initial_pressure_vt (massic_volume, temperature)               \
        , info                                                                   \
    );                                                                           \
}                                                                                \
    constexpr auto                                                               \
NAME##_vt (                                                                      \
      auto const& massic_volume                                                  \
    , auto const& temperature                                                    \
){                                                                               \
    return NAME##_vt (massic_volume, temperature, info::none);                   \
}                                                                                \
    template <info_tag_t InfoTag = info::tag::none>                              \
    constexpr auto                                                               \
NAME##_dt (                                                                      \
      auto const& density                                                        \
    , auto const& temperature                                                    \
    , auto const& initial_guess                                                  \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_vt (1. / density, temperature, initial_guess, info);           \
}                                                                                \
    template <info_tag_t InfoTag>                                                \
    constexpr auto                                                               \
NAME##_dt (                                                                      \
      auto const& density                                                        \
    , auto const& temperature                                                    \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_vt (1. / density, temperature, info);                          \
}                                                                                \
    constexpr auto                                                               \
NAME##_dt (                                                                      \
      auto const& density                                                        \
    , auto const& temperature                                                    \
){                                                                               \
    return NAME##_dt (density, temperature, info::none);                         \
}                                                                                \
    template <info_tag_t InfoTag = info::tag::none>                              \
    constexpr auto                                                               \
NAME##_tv (                                                                      \
      auto const& temperature                                                    \
    , auto const& massic_volume                                                  \
    , auto const& initial_guess                                                  \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_vt (massic_volume, temperature, initial_guess, info);          \
}                                                                                \
    template <info_tag_t InfoTag>                                                \
    constexpr auto                                                               \
NAME##_tv (                                                                      \
      auto const& temperature                                                    \
    , auto const& massic_volume                                                  \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_vt (massic_volume, temperature, info);                         \
}                                                                                \
    constexpr auto                                                               \
NAME##_tv (                                                                      \
      auto const& temperature                                                    \
    , auto const& massic_volume                                                  \
){                                                                               \
    return pressure_vt (massic_volume, temperature);                             \
}                                                                                \
    template <info_tag_t InfoTag = info::tag::none>                              \
    constexpr auto                                                               \
NAME##_td (                                                                      \
      auto const& temperature                                                    \
    , auto const& density                                                        \
    , auto const& initial_guess                                                  \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_dt (density, temperature, initial_guess, info);                \
}                                                                                \
    template <info_tag_t InfoTag>                                                \
    constexpr auto                                                               \
NAME##_td (                                                                      \
      auto const& temperature                                                    \
    , auto const& density                                                        \
    , info_t <InfoTag> info = info::none                                         \
){                                                                               \
    return NAME##_dt (density, temperature, info);                               \
}                                                                                \
    constexpr auto                                                               \
NAME##_td (                                                                      \
      auto const& temperature                                                    \
    , auto const& density                                                        \
){                                                                               \
    return NAME##_dt (density, temperature);                                     \
}                                                                                \

ISTO_IAPWS_G12_INVERSE_GEN(massic_entropy)
ISTO_IAPWS_G12_INVERSE_GEN(isothermal_compressibility)
ISTO_IAPWS_G12_INVERSE_GEN(thermal_expansion_coefficient)
ISTO_IAPWS_G12_INVERSE_GEN(massic_isobaric_heat_capacity)
ISTO_IAPWS_G12_INVERSE_GEN(massic_isochoric_heat_capacity)
ISTO_IAPWS_G12_INVERSE_GEN(speed_of_sound)
} // inline namespace g12_15
} // namespace calculisto::iapws::g12_inverse
