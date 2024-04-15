#pragma once
#include "detail/common.hpp"
#include <calculisto/root_finding/root_finding.hpp>
    using namespace calculisto::root_finding;

    namespace 
calculisto::iapws::g12
{
    inline namespace 
g12_15
{
    namespace
detail
{
        constexpr auto
    omega_0 = 0.5212269;
        constexpr auto
    L_0 = 0.76317954;
        constexpr auto
    k_0 = 0.072158686;
        constexpr auto
    k_1 = -0.31569232;
        constexpr auto
    k_2 = 5.2992608;
        constexpr auto
    T_LL = 228.2;
        constexpr auto
    rho_0 = 1081.6482;
        constexpr auto
    R = 461.523087;
        constexpr auto
    pi_0 = 300e6 / rho_0 / R / T_LL;

        constexpr auto
    c = array_t
    {
          -8.1570681381655
        ,  1.2875032
        ,  7.0901673598012
        , -3.2779161e-2
        ,  7.3703949e-1
        , -2.1628622e-1
        , -5.1782479
        ,  4.2293517e-4
        ,  2.3592109e-2
        ,  4.3773754
        , -2.9967770e-3
        , -9.6558018e-1
        ,  3.7595286
        ,  1.2632441
        ,  2.8542697e-1
        , -8.5994947e-1
        , -3.2916153e-1
        ,  9.0019616e-2
        ,  8.1149726e-2
        , -3.2788213
    };

        constexpr auto
    a = array_t
    {
           0.
        ,  0.
        ,  1.
        , -0.2555
        ,  1.5762
        ,  1.6400
        ,  3.6385
        , -0.3828
        ,  1.6219
        ,  4.3287
        ,  3.4763
        ,  5.1556
        , -0.3593
        ,  5.0361
        ,  2.9786
        ,  6.2373
        ,  4.0460
        ,  5.3558
        ,  9.0157
        ,  1.2194
    };

        constexpr auto
    b = array_t
    {
           0.
        ,  1.
        ,  0.
        ,  2.1051
        ,  1.1422
        ,  0.9510
        ,  0.
        ,  3.6402
        ,  2.0760
        , -0.0016
        ,  2.2769
        ,  0.0008
        ,  0.3706
        , -0.3975
        ,  2.9730
        , -0.3180
        ,  2.9805
        ,  2.9265
        ,  0.4456
        ,  0.1298
    };

        constexpr auto
    d = array_t
    {
           0.
        ,  0.
        ,  0.
        , -0.0016
        ,  0.6894
        ,  0.0130
        ,  0.0002
        ,  0.0435
        ,  0.0500
        ,  0.0004
        ,  0.0528
        ,  0.0147
        ,  0.8584
        ,  0.9924
        ,  1.0041
        ,  1.0961
        ,  1.0228
        ,  1.0303
        ,  1.6180
        ,  0.5213
    };

        constexpr auto
    omega (auto const& pi)
    {
        return 2. + omega_0 * pi;
    }
        using std::sqrt;
        constexpr auto
    K_2 = sqrt (1. + k_2 * k_2);

        constexpr auto
    K_1 (auto const& tau, auto const& pi)
    {
            using std::sqrt;
        return sqrt (
              pow (1. + k_0 * k_2 + k_1 * (pi - k_2 * tau), 2.) 
            - 4. * k_0 * k_1 * k_2 * (pi - k_2 * tau)
        );
    }

        constexpr auto
    L (auto const& tau, auto const& pi)
    {
        return L_0 * K_2 / 2. / k_1 / k_2 * (
              1. + k_0 * k_2
            + k_1 * (pi + k_2 * tau) - K_1 (tau, pi)
        );
    }

        constexpr auto
    L_t (auto const& tau, auto const& pi)
    {
        return (L_0 * K_2) / 2. * (1. 
            + (1. - k_0 * k_2 + k_1 * (pi - k_2 * tau)) / K_1 (tau, pi)
        );
    }

        constexpr auto
    L_p (auto const& tau, auto const& pi)
    {
        return (L_0 * K_2 
            * (K_1 (tau, pi) + k_0 * k_2 - k_1 * pi + k_1 * k_2 * tau - 1.)
        ) / (2. * k_2 * K_1 (tau, pi));
    }

        constexpr auto
    L_tt (auto const& tau, auto const& pi)
    {
        return -2. * L_0 * K_2 * k_0 * k_1 * k_2 * k_2 / pow (K_1 (tau, pi), 3.);
    }

        constexpr auto
    L_tp (auto const& tau, auto const& pi)
    {
        return 2. * L_0 * K_2 * k_0 * k_1 * k_2 / pow (K_1 (tau, pi), 3.);
    }

        constexpr auto
    L_pp (auto const& tau, auto const& pi)
    {
        return -2. * L_0 * K_2 * k_0 * k_1 / pow (K_1 (tau, pi), 3.);
    }

        constexpr auto
    Psi_r_t (auto const& tau, auto const& pi)
    {
            const auto
        tau_ = tau + 1.;
            const auto
        pi_ = pi + pi_0;
        return sum (c * a * pow (tau_, a - 1.) * pow (pi_, b) * exp (-d * pi_));
    }

        constexpr auto
    Psi_r_p (auto const& tau, auto const& pi)
    {
            const auto
        tau_ = tau + 1.;
            const auto
        pi_ = pi + pi_0;
        return sum (c * pow (tau_, a) * pow (pi_, b - 1.) * (b - d * pi_) * exp (-d * pi_));
    }

        constexpr auto
    Psi_r_tt (auto const& tau, auto const& pi)
    {
            const auto
        tau_ = tau + 1.;
            const auto
        pi_ = pi + pi_0;
        return sum (c * a * (a - 1.) * pow (tau_, a - 2.) * pow (pi_, b) * exp (-d * pi_));
    }
        constexpr auto
    Psi_r_tp (auto const& tau, auto const& pi)
    {
            const auto
        tau_ = tau + 1.;
            const auto
        pi_ = pi + pi_0;
        return sum (c * a * pow (tau_, a - 1.) * pow (pi_, b - 1.) * (b - d * pi_) * exp (-d * pi_));
    }

        constexpr auto
    Psi_r_pp (auto const& tau, auto const& pi)
    {
            const auto
        tau_ = tau + 1.;
            const auto
        pi_ = pi + pi_0;
        return sum (c * pow (tau_, a) * pow (pi_, b - 2.) * (pow (d * pi_ - b, 2.) - b) * exp (-d * pi_));
    }


        template <class T>
        constexpr auto
    get_x (T const& tau, T const& pi)
    {
            auto
        f = [&](auto const x)
        {
                using std::log;
            return 
                  L (tau, pi) 
                + log (x / (1. - x)) 
                + omega (pi) * (1. - 2. * x)
            ;
        };
            const auto
        _L = L (tau, pi);
            const auto
        _omega = omega (pi);
            T
        low;
            T
        hig;
            using std::log, std::exp, std::min;
            const auto
        a = (10. / 9.) * (log (19.) - _L);
            const auto
        b = (50. / 49.) * (log (99.) - _L);
        if (_omega < a)
        {
            low = 0.049;
            hig = 0.5;
        }
        else if (a <= _omega && _omega < b)
        {
            low = 0.0099;
            hig = 0.051;
        }
        else
        {
            low = 0.99 * exp (-(50. / 49.) * _L - _omega);
            hig = fmin (1.1 * exp (-_L - _omega), 0.0101);
        }
        return zhang (f, low, hig, { .converged = [](T a, T b, T fa, T fb)
        {
            return fa == 0 || fb == 0 || fabs (a - b) < 1e-8;
        }});
    }

        constexpr auto
    p_H (auto const& temperature)
    {
            const auto
        theta = temperature / 235.15;
        return (
              0.1 
            + 228.27 * (1. - pow (theta, 6.243))
            + 15.724 * (1. - pow (theta, 79.81))
        ) * 1e6;
    }
    /*
        constexpr auto
    dp_H_dt (auto const& temperature)
    {
            const auto
        theta = temperature / 235.15;
        return (
            + 228.27 * 6.243 * pow (theta, 6.243 - 1.)
            + 15.724 * 79.81 * pow (theta, 79.81 - 1.)
        ) / 235.15 
#ifdef ISTO_IAPWS_FLAVOR_
         / 1.
#endif
        ;
    }
    */
        constexpr auto
    T_H (auto const& pressure)
    {
            const auto
        p = pressure / (1e6);
        return (
              172.82
            + 0.03718 * p
            + 3.403e-5 * p * p
            - 1.573e-8 * p * p * p
            );
    }
} // namespace detail
    constexpr auto
massic_volume_tp (
      auto const& temperature
    , auto const& pressure
){
        using namespace detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
    return (
          (tau + 1) / 2. 
          * (omega_0 / 2. * (1 - phi * phi) + L_p (tau, pi) * (phi + 1.))
        + Psi_r_p (tau, pi)
    ) / rho_0;
}
    constexpr auto
massic_volume_pt (
      auto const& pressure
    , auto const& temperature
){
    return massic_volume_tp (temperature, pressure);
}

    constexpr auto
density_tp (
      auto const& temperature
    , auto const& pressure
){
    return 1. / massic_volume_tp (temperature, pressure);
}
    constexpr auto
density_pt (
      auto const& pressure
    , auto const& temperature
){
    return density_tp (temperature, pressure);
}

    constexpr auto
massic_entropy_tp (
      auto const& temperature
    , auto const& pressure
){
        using namespace detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
    return - R * ((tau + 1.) * L_t (tau, pi) / 2. * (phi + 1.) + (
          x * L (tau, pi) 
        + x * log (x) 
        + (1. - x) * log (1. - x) 
        + omega (pi) * x * (1. - x)
    ) + Psi_r_t (tau, pi));
}
    constexpr auto
massic_entropy_pt (
      auto const& pressure
    , auto const& temperature
){
    return massic_entropy_tp (temperature, pressure);
}

    constexpr auto
isothermal_compressibility_tp (
      auto const& temperature
    , auto const& pressure
){
        using namespace detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
        const auto
    chi = 1. / (2. / (1. - phi * phi) - omega (pi));
    return density_tp (temperature, pressure) / rho_0 / rho_0 / R / T_LL * (
        (tau + 1.) / 2. * (
              chi * pow (L_p (tau, pi) - omega_0 * phi, 2.) 
              - (phi + 1.) * L_pp (tau, pi)
        ) 
    - Psi_r_pp (tau, pi));
}
    constexpr auto
isothermal_compressibility_pt (
      auto const& pressure
    , auto const& temperature
){
    return isothermal_compressibility_tp (temperature, pressure);
}
    constexpr auto
thermal_expansion_coefficient_tp (
      auto const& temperature
    , auto const& pressure
){
        using namespace detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
        const auto
    chi = 1. / (2. / (1. - phi * phi) - omega (pi));
    return density_tp (temperature, pressure) / rho_0 / T_LL * (
          L_tp (tau, pi) / 2. * (tau + 1.) * (phi + 1.)
        + 0.5 * (omega_0 * (1. - phi * phi) / 2. + L_p (tau, pi) * (phi + 1.))
        - (tau + 1.) * L_t (tau, pi) / 2. * chi * (L_p (tau, pi) - omega_0 * phi)
        + Psi_r_tp (tau, pi)
    );
}
    constexpr auto
thermal_expansion_coefficient_pt (
      auto const& pressure
    , auto const& temperature
){
    return thermal_expansion_coefficient_tp (temperature, pressure);
}

    constexpr auto
massic_isobaric_heat_capacity_tp (
      auto const& temperature
    , auto const& pressure
){
        using namespace detail;
        const auto
    tau = temperature / T_LL - 1.;
        const auto
    pi = pressure / rho_0 / R / T_LL;
        const auto
    x = get_x (tau, pi);
        const auto
    phi = 2. * x - 1.;
        const auto
    chi = 1. / (2. / (1. - phi * phi) - omega (pi));
    return -R * (tau + 1.) * (
          L_t (tau, pi) * (phi + 1)
        + 0.5 * (tau + 1.) * (L_tt (tau, pi) * (phi + 1) - pow (L_t (tau, pi), 2.) * chi)
        + Psi_r_tt (tau, pi)
    );
}
    constexpr auto
massic_isobaric_heat_capacity_pt (
      auto const& pressure
    , auto const& temperature
){
    return massic_isobaric_heat_capacity_tp (temperature, pressure);
}
    constexpr auto
massic_isochoric_heat_capacity_tp (
      auto const& temperature
    , auto const& pressure
){
    return massic_isobaric_heat_capacity_tp (temperature, pressure)
        - temperature
        * pow (thermal_expansion_coefficient_tp (temperature, pressure), 2) 
        / density_tp (temperature, pressure) 
        / isothermal_compressibility_tp (temperature, pressure)
    ;
}
    constexpr auto
massic_isochoric_heat_capacity_pt (
      auto const& pressure
    , auto const& temperature
){
    return massic_isochoric_heat_capacity_tp (temperature, pressure);
}

    constexpr auto
speed_of_sound_tp (
      auto const& temperature
    , auto const& pressure
){
        using std::sqrt;
    return pow (
          density_tp  (temperature, pressure) 
        * isothermal_compressibility_tp (temperature, pressure) 
        * massic_isochoric_heat_capacity_tp (temperature, pressure) 
        / massic_isobaric_heat_capacity_tp (temperature, pressure)
        , - 0.5
    );
}
    constexpr auto
speed_of_sound_pt (
      auto const& pressure
    , auto const& temperature
){
    return speed_of_sound_tp (temperature, pressure);
}

    constexpr auto
homogeneous_ice_nucleation_limit_temperature_high_p (auto const& pressure)
{
    return detail::T_H (pressure);
}
    constexpr auto
homogeneous_ice_nucleation_limit_temperature_low_t (auto const& temperature)
{
    return detail::p_H (temperature);
}

    constexpr auto
homogeneous_ice_nucleation_limit_temperature_p (auto const& pressure)
{
        using namespace detail;
        const auto
    P_ = 198871388.776809;
    if (pressure > P_)
    {
        return T_H (pressure);
    }
    return zhang (
          [=](auto t){ return p_H (t) - pressure; }
        , 180.
        , 240.
        , { .converged = [](auto a, auto b, auto fa, auto fb)
            {
                return fa == 0 || fb == 0 || fabs (a - b) < 1e-10;
            }}
    );
    /*
    return newton (
          [=](auto t){ return p_H (t) * 1e6 - pressure; }
        , [=](auto t){ return dp_H_dt (t); }
        , 181.4
        , [](auto x){ return fabs (x) < 1e-8; }
    );
    */
}

} // namespace g12_15
} // namespace calculisto::iapws::g12
