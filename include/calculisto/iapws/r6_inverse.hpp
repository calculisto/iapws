#pragma once
#include "r6.hpp"
#include "r7.hpp"
#include <calculisto/root_finding/root_finding.hpp>
    using namespace calculisto::root_finding;

    namespace
calculisto::thermodynamics::iapws::r6_inverse
{
    inline namespace
r6_95_2016
{
    constexpr auto
initial_density_pt (
      auto const& pressure
    , auto const& temperature
){
    return (
            pressure > 100e6 
        || (
               temperature > 1073.15 
            && pressure > 50e6
        ))
    ? 1000.
    : r7::density_pt (pressure, temperature)
;
}

    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
density_pt (
      auto const& pressure
    , auto const& temperature
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
        using namespace r6;
        using namespace r6::detail;
        auto
    tau = critical_temperature / temperature;
    return newton (
          [=](auto density)
          { 
            return r6::pressure_dt (density, temperature) - pressure; 
          }
        , [=](auto density)
          {
                using namespace detail;
                auto
            delta = density / critical_density;
            return (1 + 2 * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau)) * massic_gas_constant * temperature;
          }
        , initial_guess
        , { 
            .converged = [](auto curr, auto prev, auto f)
            { 
                return fabs (f) < 1e-8 || fabs ((curr - prev)/curr) < 1e-8; 
            } 
          } // options
        , info
    );
}
    template <info_tag_t InfoTag>
    constexpr auto
density_pt (
      auto const& pressure
    , auto const& temperature
    , info_t <InfoTag> info
){
    return density_pt (pressure, temperature, initial_density_pt (pressure, temperature), info);
}
    constexpr auto
density_pt (
      auto const& pressure
    , auto const& temperature
){
    return density_pt (pressure, temperature, initial_density_pt (pressure, temperature), info::none);
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
density_tp (
      auto const& temperature
    , auto const& pressure
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return density_pt (pressure, temperature, initial_guess, info);
}
    template <info_tag_t InfoTag>
    constexpr auto
density_tp (
      auto const& temperature
    , auto const& pressure
    , info_t <InfoTag> info
){
    return density_pt (pressure, temperature, info);
}
    constexpr auto
density_tp (
      auto const& temperature
    , auto const& pressure
){
    return density_pt (pressure, temperature, info::none);
}

#define ISTO_IAPWS_R6_INVERSE_GEN(NAME)                                      \
    template <info_tag_t InfoTag = info::tag::none>                          \
    constexpr auto                                                           \
NAME##_pt (                                                                  \
      auto const& pressure                                                   \
    , auto const& temperature                                                \
    , auto const& density_initial_guess                                      \
    , info_t <InfoTag> info = info::none                                     \
){                                                                           \
        const auto                                                           \
    d = density_pt (pressure, temperature, density_initial_guess, info);     \
    if constexpr (InfoTag != info::tag::none)                                \
    {                                                                        \
        return std::pair { r6::NAME##_dt (d.first, temperature), d.second }; \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        return r6::NAME##_dt (d, temperature);                               \
    }                                                                        \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag>                                            \
    constexpr auto                                                           \
NAME##_pt (                                                                  \
      auto const& pressure                                                   \
    , auto const& temperature                                                \
    , info_t <InfoTag> info                                                  \
){                                                                           \
    return NAME##_pt (                                                       \
          pressure                                                           \
        , temperature                                                        \
        , initial_density_pt (pressure, temperature)                         \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    constexpr auto                                                           \
NAME##_pt (                                                                  \
      auto const& pressure                                                   \
    , auto const& temperature                                                \
){                                                                           \
    return NAME##_pt (                                                       \
          pressure                                                           \
        , temperature                                                        \
        , initial_density_pt (pressure, temperature)                         \
        , info::none                                                         \
    );                                                                       \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag = info::tag::none>                          \
    constexpr auto                                                           \
NAME##_tp (                                                                  \
      auto const& temperature                                                \
    , auto const& pressure                                                   \
    , auto const& density_initial_guess                                      \
    , info_t <InfoTag> info = info::none                                     \
){                                                                           \
    return NAME##_pt (                                                       \
          pressure                                                           \
        , temperature                                                        \
        , density_initial_guess                                              \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag>                                            \
    constexpr auto                                                           \
NAME##_tp (                                                                  \
      auto const& temperature                                                \
    , auto const& pressure                                                   \
    , info_t <InfoTag> info                                                  \
){                                                                           \
    return NAME##_pt (                                                       \
          pressure                                                           \
        , temperature                                                        \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    constexpr auto                                                           \
NAME##_tp (                                                                  \
      auto const& temperature                                                \
    , auto const& pressure                                                   \
){                                                                           \
    return NAME##_pt (                                                       \
          pressure                                                           \
        , temperature                                                        \
        , info::none                                                         \
    );                                                                       \
}                                                                                      

ISTO_IAPWS_R6_INVERSE_GEN(massic_internal_energy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_entropy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_enthalpy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_isochoric_heat_capacity)
ISTO_IAPWS_R6_INVERSE_GEN(massic_isobaric_heat_capacity)
//ISTO_IAPWS_R6_INVERSE_GEN(joule_thompson_coefficient)
ISTO_IAPWS_R6_INVERSE_GEN(massic_gibbs_free_energy)
ISTO_IAPWS_R6_INVERSE_GEN(speed_of_sound)
ISTO_IAPWS_R6_INVERSE_GEN(isothermal_stress_coefficient)
ISTO_IAPWS_R6_INVERSE_GEN(relative_pressure_coefficient)
//ISTO_IAPWS_R6_INVERSE_GEN(isobaric_cubic_expansion_coefficient)
//ISTO_IAPWS_R6_INVERSE_GEN(isothermal_compressibility)
#undef ISTO_IAPWS_R6_INVERSE_GEN

// Saturation line
    namespace
detail
{
        template <class T>
        auto
    f_saturation (Eigen::Matrix <T, 3, 1> const& x, auto const tau)
    {
            using namespace r6::detail;
            auto
        r = Eigen::Matrix <T, 3, 1> {};
            const auto&
        pi = x[0];
            const auto&
        delta_p = x[1];
            const auto&
        delta_pp = x[2];
        r[0] = delta_p * (1 + delta_p * phi_r_d (delta_p, tau)) - pi;
        r[1] = delta_pp * (1 + delta_pp * phi_r_d (delta_pp, tau)) - pi;
        r[2] = 
              phi_r (delta_p, tau) 
            - phi_r (delta_pp, tau) 
            - pi * (1 / delta_pp - 1 / delta_p)
            + log (delta_p / delta_pp)
        ;
        return r;
    }
} // namespace detail

    template <
          class T
        , info_tag_t InfoTag
        , class Range = Eigen::Matrix <T, 3, 1>
    >
    auto
saturation_pressure_t (
      T const& temperature
    , T const& pressure_initial_guess
    , T const& density_liquid_initial_guess
    , T const& density_gas_initial_guess
    , multidimensional_newton_options_t <3, Range, Range, Range> const& options
    , info_t <InfoTag> info
){
        const auto
    pi_0 = pressure_initial_guess / massic_gas_constant / temperature / critical_density;
        const auto
    delta_p_0 = density_liquid_initial_guess / critical_density;
        const auto
    delta_pp_0 = density_gas_initial_guess / critical_density;
        const auto
    x_0 = Eigen::Matrix <T, 3, 1> {
          pi_0
        , delta_p_0
        , delta_pp_0
    };
        const auto
    tau = critical_temperature / temperature;
        const auto
    f = [tau](auto const& x){ return detail::f_saturation (x, tau); };
    if constexpr (InfoTag == info::tag::none)
    {
            const auto
        r = newton <3> (
              f
            , x_0
            , options
            , info
        );
        return std::tuple {
              r[0] * massic_gas_constant * temperature * critical_density
            , r[1] * critical_density
            , r[2] * critical_density
        };
    }
    else
    {
            const auto
        [r, info_data] = newton <3> (
              f
            , x_0
            , options
            , info
        );
        return std::tuple {
              r[0] * massic_gas_constant * temperature * critical_density
            , r[1] * critical_density
            , r[2] * critical_density
            , info_data
        };
    }
}
    template <
          class T
        , info_tag_t InfoTag
        , class Range = Eigen::Matrix <T, 3, 1>
    >
    auto
saturation_pressure_t (
      T const& temperature
    , multidimensional_newton_options_t <3, Range, Range, Range> const& options
    , info_t <InfoTag> info
){
        const auto
    pressure_inital_guess = r7::saturation_pressure_t (temperature);
    return saturation_pressure_t (
          temperature
        , pressure_inital_guess
        , r7::r1::density_pt (pressure_inital_guess, temperature)
        , r7::r2::density_pt (pressure_inital_guess, temperature)
        , options
        , info
    );
}
    template <
          class T
        , info_tag_t InfoTag = info::tag::none
        , class Range = Eigen::Matrix <T, 3, 1>
    >
    auto
saturation_pressure_t (
      T const& temperature
    , T const& relative_tolerance = 1e-8
    , T const& absolute_tolerance = 0.
    , int max_iter = 100
    , info_t <InfoTag> info = info::none
){
        const auto
    pressure_inital_guess = r7::saturation_pressure_t (temperature);
    return saturation_pressure_t (
          temperature
        , pressure_inital_guess
        , r7::r1::density_pt (pressure_inital_guess, temperature)
        , r7::r2::density_pt (pressure_inital_guess, temperature)
        , {
              .max_iter = max_iter
            , .converged = [&](
                  Eigen::Matrix <T, 3, 1> const& current
                , Eigen::Matrix <T, 3, 1> const& past
                , Eigen::Matrix <T, 3, 1> const& result
              ){
                return
                       (past - current).norm () / current.norm ()
                       <
                       relative_tolerance
                    || result.squaredNorm () == absolute_tolerance
                ;
              }
          }
        , info
    );
}

    constexpr auto
saturation_temperature_p (auto const& pressure)
{
    // FIXME:
    throw not_yet_implemented_e {};
}

// Temperature
    constexpr auto
initial_temperature_pd (
      auto const& pressure [[maybe_unused]]
    , auto const& density  [[maybe_unused]]
){
    return critical_temperature;
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
temperature_pd (
      auto const& pressure
    , auto const& density
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
        using namespace r6;
        using namespace r6::detail;
        auto
    delta = density / critical_density;
    return newton (
          [=](auto temperature)
          { 
            return r6::pressure_dt (density, temperature) - pressure; 
          }
        , [=](auto temperature)
          {
                using namespace detail;
                auto
            tau = critical_temperature / temperature;
            return -tau * density * massic_gas_constant * (-1 / tau - delta / tau * phi_r_d (delta, tau) + delta * phi_r_dt (delta, tau));
          }
        , initial_guess
        , { 
            .converged = [](auto curr, auto prev, auto f)
            { 
                return fabs (f) < 1e-8 || fabs ((curr - prev)/curr) < 1e-8; 
            } 
          } // options
        , info
    );
}
    template <info_tag_t InfoTag>
    constexpr auto
temperature_pd (
      auto const& pressure
    , auto const& density
    , info_t <InfoTag> info
){
    return temperature_pd (pressure, density, initial_temperature_pd (pressure, density), info);
}
    constexpr auto
temperature_pd (
      auto const& pressure
    , auto const& density
){
    return temperature_pd (pressure, density, initial_temperature_pd (pressure, density), info::none);
}
    template <info_tag_t InfoTag = info::tag::none>
    constexpr auto
temperature_dp (
      auto const& density
    , auto const& pressure
    , auto const& initial_guess
    , info_t <InfoTag> info = info::none
){
    return temperature_pd (pressure, density, initial_guess, info);
}
    template <info_tag_t InfoTag>
    constexpr auto
temperature_dp (
      auto const& density
    , auto const& pressure
    , info_t <InfoTag> info
){
    return temperature_pd (pressure, density, info);
}
    constexpr auto
temperature_dp (
      auto const& density
    , auto const& pressure
){
    return temperature_pd (pressure, density, info::none);
}

#define ISTO_IAPWS_R6_INVERSE_GEN(NAME)                                      \
    template <info_tag_t InfoTag = info::tag::none>                          \
    constexpr auto                                                           \
NAME##_pd (                                                                  \
      auto const& pressure                                                   \
    , auto const& density                                                    \
    , auto const& temperature_initial_guess                                  \
    , info_t <InfoTag> info = info::none                                     \
){                                                                           \
        const auto                                                           \
    t = temperature_pd (pressure, density, temperature_initial_guess, info); \
    if constexpr (InfoTag != info::tag::none)                                \
    {                                                                        \
        return std::pair { r6::NAME##_dt (density, t.first), t.second };     \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        return r6::NAME##_dt (density, t);                                   \
    }                                                                        \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag>                                            \
    constexpr auto                                                           \
NAME##_pd (                                                                  \
      auto const& pressure                                                   \
    , auto const& density                                                    \
    , info_t <InfoTag> info                                                  \
){                                                                           \
    return NAME##_pd (                                                       \
          pressure                                                           \
        , density                                                            \
        , initial_temperature (pressure, density)                            \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    constexpr auto                                                           \
NAME##_pd (                                                                  \
      auto const& pressure                                                   \
    , auto const& density                                                    \
){                                                                           \
    return NAME##_pd (                                                       \
          pressure                                                           \
        , density                                                            \
        , initial_temperature (pressure, density)                            \
        , info::none                                                         \
    );                                                                       \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag = info::tag::none>                          \
    constexpr auto                                                           \
NAME##_dp (                                                                  \
      auto const& density                                                    \
    , auto const& pressure                                                   \
    , auto const& temperature_initial_guess                                  \
    , info_t <InfoTag> info = info::none                                     \
){                                                                           \
    return NAME##_pd (                                                       \
          pressure                                                           \
        , density                                                            \
        , temperature_initial_guess                                          \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    template <info_tag_t InfoTag>                                            \
    constexpr auto                                                           \
NAME##_dp (                                                                  \
      auto const& density                                                    \
    , auto const& pressure                                                   \
    , info_t <InfoTag> info                                                  \
){                                                                           \
    return NAME##_pd (                                                       \
          pressure                                                           \
        , density                                                            \
        , info                                                               \
    );                                                                       \
}                                                                            \
                                                                             \
    constexpr auto                                                           \
NAME##_dp (                                                                  \
      auto const& density                                                    \
    , auto const& pressure                                                   \
){                                                                           \
    return NAME##_pd (                                                       \
          pressure                                                           \
        , density                                                            \
        , info::none                                                         \
    );                                                                       \
}                                                                                      

ISTO_IAPWS_R6_INVERSE_GEN(massic_internal_energy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_entropy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_enthalpy)
ISTO_IAPWS_R6_INVERSE_GEN(massic_isochoric_heat_capacity)
ISTO_IAPWS_R6_INVERSE_GEN(massic_isobaric_heat_capacity)
//ISTO_IAPWS_R6_INVERSE_GEN(joule_thompson_coefficient)
ISTO_IAPWS_R6_INVERSE_GEN(massic_gibbs_free_energy)
ISTO_IAPWS_R6_INVERSE_GEN(speed_of_sound)
ISTO_IAPWS_R6_INVERSE_GEN(isothermal_stress_coefficient)
ISTO_IAPWS_R6_INVERSE_GEN(relative_pressure_coefficient)
//ISTO_IAPWS_R6_INVERSE_GEN(isobaric_cubic_expansion_coefficient)
//ISTO_IAPWS_R6_INVERSE_GEN(isothermal_compressibility)
#undef ISTO_IAPWS_R6_INVERSE_GEN

} // inline namespace r6_95_2016
} // namespace calculisto::thermodynamics::iapws::r6_inverse
