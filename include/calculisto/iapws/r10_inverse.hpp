#pragma once
#include "r10.hpp"
#include <calculisto/root_finding/root_finding.hpp>
    using namespace calculisto::root_finding;

    namespace
calculisto::iapws::r10_inverse
{
    namespace
detail
{
    template <class T>
    constexpr auto
initial_pressure_vt (
      [[maybe_unused]] T const& massic_volume
    , [[maybe_unused]] T const& temperature
){
    return T { 1e5 };
}
} // namespace detail
    inline namespace
r10_06_2009
{
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
pressure_vt (
      auto const& massic_volume
    , auto const& temperature
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
        using namespace r10::detail;
    return newton (
          [=](auto pressure)
          { 
                return g_p (pressure, temperature) - massic_volume; 
          }
        , [=](auto pressure)
          {
            return g_pp (pressure, temperature) ;
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


#define ISTO_IAPWS_R10_INVERSE_GEN(NAME)                                         \
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
        return std::pair { r10::NAME##_pt (p.first, temperature), p.second };    \
    }                                                                            \
    else                                                                         \
    {                                                                            \
        return r10::NAME##_pt (p, temperature);                                  \
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

ISTO_IAPWS_R10_INVERSE_GEN(massic_entropy)
ISTO_IAPWS_R10_INVERSE_GEN(massic_isobaric_heat_capacity)
ISTO_IAPWS_R10_INVERSE_GEN(massic_enthalpy)
ISTO_IAPWS_R10_INVERSE_GEN(massic_internal_energy)
ISTO_IAPWS_R10_INVERSE_GEN(massic_helmholtz_energy)
ISTO_IAPWS_R10_INVERSE_GEN(cubic_expansion_coefficient)
ISTO_IAPWS_R10_INVERSE_GEN(pressure_coefficient)
ISTO_IAPWS_R10_INVERSE_GEN(isothermal_compressibility)
ISTO_IAPWS_R10_INVERSE_GEN(isentropic_compressibility)
} // inline namespace r10_06_2009
} // namespace calculisto::iapws::r10_inverse
