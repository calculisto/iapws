#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r6.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::iapws;
    using namespace calculisto::iapws::r6;
#include "../include/calculisto/iapws/detail/data_for_the_tests.hpp"

TEST_CASE("r6.hpp")
{
SUBCASE("base functions")
{
    CHECK(detail::phi_0    (838.025 / 322., 647.096 / 500.) == Approx {  0.204797734e1 }.epsilon (1e-8));
    CHECK(detail::phi_r    (838.025 / 322., 647.096 / 500.) == Approx { -0.342693206e1 }.epsilon (1e-8));
    CHECK(detail::phi_0_d  (838.025 / 322., 647.096 / 500.) == Approx {  0.384236747   }.epsilon (1e-8));
    CHECK(detail::phi_r_d  (838.025 / 322., 647.096 / 500.) == Approx { -0.364366650   }.epsilon (1e-8));
    CHECK(detail::phi_0_dd (838.025 / 322., 647.096 / 500.) == Approx { -0.147637878   }.epsilon (1e-8));
    CHECK(detail::phi_r_dd (838.025 / 322., 647.096 / 500.) == Approx {  0.856063701   }.epsilon (1e-8));
    CHECK(detail::phi_0_t  (838.025 / 322., 647.096 / 500.) == Approx {  0.904611106e1 }.epsilon (1e-8));
    CHECK(detail::phi_r_t  (838.025 / 322., 647.096 / 500.) == Approx { -0.581403435e1 }.epsilon (1e-8));
    CHECK(detail::phi_0_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.193249185e1 }.epsilon (1e-8));
    CHECK(detail::phi_r_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.223440737e1 }.epsilon (1e-8));
    //CHECK(detail::phi_0_dt (838.025 / 322., 647.096 / 500.) == Approx {  0.            }.epsilon (1e-8));
    CHECK(detail::phi_r_dt (838.025 / 322., 647.096 / 500.) == Approx { -0.112176915e1 }.epsilon (1e-8));
    CHECK(detail::phi_0    (358.000 / 322., 647.096 / 647.) == Approx { -0.156319605e1 }.epsilon (1e-8));
    CHECK(detail::phi_r    (358.000 / 322., 647.096 / 647.) == Approx { -0.121202657e1 }.epsilon (1e-8));
    CHECK(detail::phi_0_d  (358.000 / 322., 647.096 / 647.) == Approx {  0.899441341   }.epsilon (1e-8));
    CHECK(detail::phi_r_d  (358.000 / 322., 647.096 / 647.) == Approx { -0.714012024   }.epsilon (1e-8));
    CHECK(detail::phi_0_dd (358.000 / 322., 647.096 / 647.) == Approx { -0.808994726   }.epsilon (1e-8));
    CHECK(detail::phi_r_dd (358.000 / 322., 647.096 / 647.) == Approx {  0.475730696   }.epsilon (1e-8));
    CHECK(detail::phi_0_t  (358.000 / 322., 647.096 / 647.) == Approx {  0.980343918e1 }.epsilon (1e-8));
    CHECK(detail::phi_r_t  (358.000 / 322., 647.096 / 647.) == Approx { -0.321722501e1 }.epsilon (1e-8));
    CHECK(detail::phi_0_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.343316334e1 }.epsilon (1e-8));
    CHECK(detail::phi_r_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.996029507e1 }.epsilon (1e-8));
    //CHECK(detail::phi_0_dt (358.000 / 322., 647.096 / 647.) == Approx {  0.            }.epsilon (1e-8));
    CHECK(detail::phi_r_dt (358.000 / 322., 647.096 / 647.) == Approx { -0.133214720e1 }.epsilon (1e-8));

} // SUBCASE("base functions")
SUBCASE("Third-order derivatives")
{
    // The analytical expressions were established by us. We test them against
    // finite difference.
        auto
    fd_t = [](auto f, double delta, double tau, double d = 1e-6)
    {
        return (
            - 0.5 * f (delta, tau - d)
            + 0.5 * f (delta, tau + d)
        ) / d;
    };
        auto
    fd_d = [](auto f, double delta, double tau, double d = 1e-6)
    {
        return (
            - 0.5 * f (delta - d, tau)
            + 0.5 * f (delta + d, tau)
        ) / d;
    };
        auto
    fd_tt = [](auto f, double delta, double tau, double d = 1e-4)
    {
        return (
                  f (delta, tau - d)
            - 2 * f (delta, tau)
            +     f (delta, tau + d)
        ) / d / d;
    };
        auto
    fd_dd = [](auto f, double delta, double tau, double d = 1e-4)
    {
        return (
                  f (delta - d, tau)
            - 2 * f (delta, tau)
            +     f (delta + d, tau)
        ) / d / d;
    };
        auto
    fd_dt = [](auto f, double delta, double tau, double d = 1e-4)
    {
        return (
              f (delta - d, tau - d)
            - f (delta - d, tau + d)
            - f (delta + d, tau - d)
            + f (delta + d, tau + d)
        ) / 4 / d / d;
    };
        auto
    fd_ddd = [](auto f, double delta, double tau, double d = 1e-3)
    {
        return (
            - 0.5 * f (delta - d - d, tau)
            +       f (delta - d, tau)
            -       f (delta + d, tau)
            + 0.5 * f (delta + d + d, tau)
        ) / d / d / d;
    };
        auto
    fd_ttt = [](auto f, double delta, double tau, double d = 1e-3)
    {
        return (
            - 0.5 * f (delta, tau - d - d)
            +       f (delta, tau - d)
            -       f (delta, tau + d)
            + 0.5 * f (delta, tau + d + d)
        ) / d / d / d;
    };
        auto
    fd_ddt = [](auto f, double delta, double tau, double d = 1e-3)
    {
        return (
            - 0.5 * f (delta - d, tau - d)
            +       f (delta, tau - d)
            - 0.5 * f (delta + d, tau - d)
            + 0.5 * f (delta - d, tau + d)
            -       f (delta, tau + d)
            + 0.5 * f (delta + d, tau + d)
        ) / d / d / d;
    };
        auto
    fd_dtt = [](auto f, double delta, double tau, double d = 1e-3)
    {
        return (
            - 0.5 * f (delta - d, tau - d)
            +       f (delta - d, tau)
            - 0.5 * f (delta - d, tau + d)
            + 0.5 * f (delta + d, tau - d)
            -       f (delta + d, tau)
            + 0.5 * f (delta + d, tau + d)
        ) / d / d / d;
    };

        const auto
    d_and_t = std::vector
    {
            std::pair
          { 838.025 / 322., 647.096 / 500. }
        , { 358.000 / 322., 647.096 / 647. }
    };

    for (auto [d, t]: d_and_t)
    {
        INFO("d= ", d, ", t= ", t);
        CHECK(detail::phi_r_d (d, t)   == Approx { fd_d   (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_d (d, t)  )).epsilon (1e-5));
        CHECK(detail::phi_r_t (d, t)   == Approx { fd_t   (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_t (d, t)  )).epsilon (1e-5));
        CHECK(detail::phi_r_dd (d, t)  == Approx { fd_dd  (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_dd (d, t) )).epsilon (1e-5));
        CHECK(detail::phi_r_tt (d, t)  == Approx { fd_tt  (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_tt (d, t) )).epsilon (1e-5));
        CHECK(detail::phi_r_ddd (d, t) == Approx { fd_ddd (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_ddd (d, t))).epsilon (1e-5));

        CHECK(detail::phi_r_dt (d, t) == Approx { fd_dt (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_dt (d, t))).epsilon (1e-5));
        CHECK(detail::phi_r_ddt (d, t) == Approx { fd_ddt (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_dt (d, t))).epsilon (1e-5));
        CHECK(detail::phi_r_dtt (d, t) == Approx { fd_dtt (&detail::phi_r <double, double>, d, t) }.scale (fabs (detail::phi_r_dt (d, t))).epsilon (1e-5));
    }

/*
    CHECK(detail::phi_r_ddd (838.025 / 322., 647.096 / 500.) == Approx { fd_ddd (&detail::phi_r <double, double>, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddd (358.000 / 322., 647.096 / 647.) == Approx { fd_ddd (&detail::phi_r <double, double>, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddt (838.025 / 322., 647.096 / 500.) == Approx { fd_ddt (&detail::phi_r <double, double>, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddt (358.000 / 322., 647.096 / 647.) == Approx { fd_ddt (&detail::phi_r <double, double>, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-5));
    CHECK(detail::phi_r_dtt (838.025 / 322., 647.096 / 500.) == Approx { fd_dtt (&detail::phi_r <double, double>, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_dtt (358.000 / 322., 647.096 / 647.) == Approx { fd_dtt (&detail::phi_r <double, double>, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-5));
    */
    /*
        auto
    d = 1e-6;
        auto
    p_ddd = [d](auto delta, auto tau)
    {
        return (
              detail::phi_r_dd (delta + d, tau)
            - detail::phi_r_dd (delta - d, tau)
        ) / 2 / d;
    };
        auto
    p_ddt = [d](auto delta, auto tau)
    {
        return (
              detail::phi_r_dd (delta, tau + d)
            - detail::phi_r_dd (delta, tau - d)
        ) / 2 / d;
    };
        auto
    p_dtt = [d](auto delta, auto tau)
    {
        return (
              detail::phi_r_dt (delta, tau + d)
            - detail::phi_r_dt (delta, tau - d)
        ) / 2 / d;
    };
    CHECK(detail::phi_r_ddd (838.025 / 322., 647.096 / 500.) == Approx { p_ddd (838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddd (358.000 / 322., 647.096 / 647.) == Approx { p_ddd (358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddt (838.025 / 322., 647.096 / 500.) == Approx { p_ddt (838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_ddt (358.000 / 322., 647.096 / 647.) == Approx { p_ddt (358.000 / 322., 647.096 / 647.) }.epsilon (1e-5));
    CHECK(detail::phi_r_dtt (838.025 / 322., 647.096 / 500.) == Approx { p_dtt (838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(detail::phi_r_dtt (358.000 / 322., 647.096 / 647.) == Approx { p_dtt (358.000 / 322., 647.096 / 647.) }.epsilon (1e-5));
    */
}
SUBCASE("Derivatives of the isobaric_cubic_expansion_coefficient")
{
        using namespace detail;
        /*
        auto
    d_alpha_v_d_tau = [](auto delta, auto tau)
    { return
        -delta * tau * phi_r_dtt (delta, tau) * tau / critical_temperature
            / (1 + 2 * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau))
        + isobaric_cubic_expansion_coefficient_dt (delta * critical_density, critical_temperature / tau)
        * (1 / tau - (2 * delta * phi_r_dt (delta, tau) + delta * delta *
                    * phi_r_ddt (delta, tau)) / (1 + 2 * delta * phi_r_d
           * (delta, tau) + delta * delta * phi_r_dd (delta, tau)));
    };
    */
        const auto
    d = 1e-6;
        auto
    d_a_v_d_t = [d](auto density, auto temperature)
    {
        return (
              isobaric_cubic_expansion_coefficient_dt (density, temperature + d)
            - isobaric_cubic_expansion_coefficient_dt (density, temperature - d)
        ) / 2 / d;
    };
        /*
    CHECK(detail::d_alpha_v_d_tau (838.025 / 322., 647.096 / 500.) == Approx { d_a_v_d_t (838.025, 500.) }.epsilon (1e-6));
    CHECK(detail::d_alpha_v_d_tau (358.000 / 322., 647.096 / 647.) == Approx { d_a_v_d_t (358.000, 647.) }.epsilon (1e-6));
    CHECK(detail::d_alpha_v_d_tau (838.025 / 322., 647.096 / 500.) == Approx { d_a_v_d_t (838.025, 500.) }.epsilon (1e-6));
    CHECK(detail::d_alpha_v_d_tau (358.000 / 322., 647.096 / 647.) == Approx { d_a_v_d_t (358.000, 647.) }.epsilon (1e-5));
    */
}
/*
SUBCASE("WIP")
{
        using namespace detail;
        auto
    s1_dd = [](auto delta, auto tau)
    {
        return sum (n_1 * d_1 * (d_1 - 1) * pow (delta, d_1 - 2) * pow (tau, t_1));
    };
        auto
    s2_dd = [](auto delta, auto tau)
    {
        return
        + sum (
              n_2 * exp (-pow (delta, c_2)) * (
              pow (delta, d_2 - 2) * pow (tau, t_2) * (
                  (d_2 - c_2 * pow (delta, c_2))
                * (d_2 - 1 - c_2 * pow (delta, c_2))
                - pow (c_2, 2) * pow (delta, c_2)
              ))
          )
        ;
    };
        auto
    s3_dd = [](auto delta, auto tau)
    {
        return
        + sum (
              n_3 * pow (tau, t_3)
            * exp(
                  -alpha * pow (delta - epsilon, 2)
                - beta_1 * pow (tau - detail::gamma, 2)
              ) * (
                  -2 * alpha * pow (delta, d_3)
                + 4 * pow (alpha, 2) * pow (delta, d_3)
                    * pow (delta - epsilon, 2)
                - 4 * d_3 * alpha * pow (delta, d_3 - 1) * (delta - epsilon)
                + d_3 * (d_3 - 1) * pow (delta, d_3 - 2)
              )
          )
        ;
    };
        auto
    s4_dd = [](auto delta, auto tau)
    {
        return
        + sum (
              n_4 * (pow (Delta(delta, tau), b) * (
                  2 * Psi_d (delta, tau)
                + delta * Psi_dd (delta, tau)
              )
            + 2 * Delta_b_d (delta, tau) * (
                  Psi (delta, tau)
                + delta * Psi_d (delta, tau)
              )
            + Delta_b_dd (delta, tau) * delta * Psi (delta, tau))
          )
        ;
    };
        auto
    s1_ddt = [](auto delta, auto tau)
    {
        return
          sum (
              n_1 * d_1 * (d_1 - 1) * t_1
            * pow (delta, d_1 - 2) * pow (tau, t_1 - 1)
          )
        ;
    };
        auto
    s2_ddt = [](auto delta, auto tau)
    {
        return
        + sum (
              n_2 * t_2 * exp (-pow (delta, c_2))
            * pow (tau, t_2 - 1) * pow (delta, d_2 - 2)
            * (
                  d_2 * d_2
                - d_2
                + pow (delta, c_2) * c_2 * (1 - c_2 - 2 * d_2)
                + pow (delta, 2 * c_2) * c_2 * c_2
              )
          )
        ;
    };
        auto
    s3_ddt = [](auto delta, auto tau)
    {
            const auto
        fd = f_d (delta, tau);
            const auto
        fdd = f_dd (delta, tau);
            const auto
        ft = f_t (delta, tau);
        return
        + sum (
              n_3 * pow (delta, d_3 - 2) * pow (tau, t_3 - 1)
            * exp (f (delta, tau))
            * (
                  ft * d_3 * d_3 * tau
                + 2 * delta * fd * ft * d_3 * tau
                - ft * d_3 * tau
                + delta * delta * fdd * ft * tau
                + delta * delta * fd * fd * ft * tau
                + d_3 * d_3 * t_3
                + 2 * delta * fd * d_3 * t_3
                - d_3 * t_3
                + delta * delta * fdd * t_3
                + delta * delta * fd * fd * t_3
              )
          )
        ;
    };
        auto
    s4_ddt = [](auto delta, auto tau)
    {
            const auto
        D = Delta (delta, tau);
            const auto
        Dd = Delta_d (delta, tau);
            const auto
        Ddd = Delta_dd (delta, tau);
            const auto
        Dddt = Delta_ddt (delta, tau);
            const auto
        Dt = Delta_t (delta, tau);
            const auto
        Ddt = Delta_dt (delta, tau);
            const auto
        P = Psi (delta, tau);
            const auto
        Pd = Psi_d (delta, tau);
            const auto
        Pdd = Psi_dd (delta, tau);
            const auto
        Pddt = Psi_ddt (delta, tau);
            const auto
        Pt = Psi_t (delta, tau);
            const auto
        Pdt = Psi_dt (delta, tau);
        return
        + sum (n_4 * pow (D, b - 3) * (
              Dd * Dd * Dt * P * delta * b * b * b
            + D * Dd * Dd * Pt * delta * b * b
            + 2 * D * Dd * Dt * Pd * delta * b * b
            + D * Ddd * Dt * P * delta * b * b
            - 3 * Dd * Dd * Dt * P * delta * b * b
            + 2 * D * Dd * Ddt * P * delta * b * b
            + 2 * D * Dd * Dt * P * b * b
            + D * D * Ddd * Pt * delta * b
            - D * Dd * Dd * Pt * delta * b
            + D * D * Dt * Pdd * delta * b
            + 2 * D * D * Dd * Pdt * delta * b
            - 2 * D * Dd * Dt * Pd * delta * b
            + 2 * D * D * Ddt * Pd * delta * b
            - D * Ddd * Dt * P * delta * b
            + 2 * Dd * Dd * Dt * P * delta * b
            + D * D * Dddt * P * delta * b
            - 2 * D * Dd * Ddt * P * delta * b
            + 2 * D * D * Dd * Pt * b
            + 2 * D * D * Dt * Pd * b
            - 2 * D * Dd * Dt * P * b
            + 2 * D * D * Ddt * P * b
            + D * D * D * Pddt * delta
            + 2 * D * D * D * Pdt
          ))
        ;
    };
        const auto
    d = 1e-5;
        auto
    dt = [d](auto f, auto delta, auto tau)
    {
        return (f (delta, tau + d) - f (delta, tau - d)) / 2 / d;
    };
        auto
    dd = [d](auto f, auto delta, auto tau)
    {
        return (f (delta + d, tau) - f (delta - d, tau)) / 2 / d;
    };
    CHECK(s1_ddt (838.025 / 322., 647.096 / 500.) == Approx { dt (s1_dd, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(s1_ddt (358.000 / 322., 647.096 / 647.) == Approx { dt (s1_dd, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));
    CHECK(s2_ddt (838.025 / 322., 647.096 / 500.) == Approx { dt (s2_dd, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(s2_ddt (358.000 / 322., 647.096 / 647.) == Approx { dt (s2_dd, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));
    CHECK(s3_ddt (838.025 / 322., 647.096 / 500.) == Approx { dt (s3_dd, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(s3_ddt (358.000 / 322., 647.096 / 647.) == Approx { dt (s3_dd, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));
    CHECK(s4_ddt (838.025 / 322., 647.096 / 500.) == Approx { dt (s4_dd, 838.025 / 322., 647.096 / 500.) }.epsilon (1e-6));
    CHECK(s4_ddt (358.000 / 322., 647.096 / 647.) == Approx { dt (s4_dd, 358.000 / 322., 647.096 / 647.) }.epsilon (1e-6));

#define YOP(FUN, DER) \
CHECK(sum (DER (838.025 / 322., 647.096 / 500.)) == Approx { sum (dd (FUN <double, double>, 838.025 / 322., 647.096 / 500.)) }.epsilon (1e-5)); \
CHECK(sum (DER (358.000 / 322., 647.096 / 647.)) == Approx { sum (dd (FUN <double, double>, 358.000 / 322., 647.096 / 647.)) }.epsilon (1e-5));

#define YIP(FUN, DER) \
CHECK(sum (DER (838.025 / 322., 647.096 / 500.)) == Approx { sum (dt (FUN <double, double>, 838.025 / 322., 647.096 / 500.)) }.epsilon (1e-5)); \
CHECK(sum (DER (358.000 / 322., 647.096 / 647.)) == Approx { sum (dt (FUN <double, double>, 358.000 / 322., 647.096 / 647.)) }.epsilon (1e-5));

    YOP(Theta, Theta_d)
    YOP(Theta_d, Theta_dd)
    YOP(Theta_dd, Theta_ddd)

    YOP(Psi, Psi_d)
    YOP(Psi_d, Psi_dd)
    YOP(Psi_dd, Psi_ddd)
    YIP(Psi, Psi_t)
    YIP(Psi_d, Psi_dt)
    YIP(Psi_dd, Psi_ddt)

    YOP(Delta, Delta_d)
    YOP(Delta_d, Delta_dd)
    YOP(Delta_dd, Delta_ddd)
    YIP(Delta, Delta_t)
    YIP(Delta_d, Delta_dt)
    YIP(Delta_dd, Delta_ddt)

    YOP(f, f_d)
    YOP(f_d, f_dd)
    YIP(f, f_t)
}
*/
SUBCASE("main API")
{
    for(const auto& e: calculisto::iapws::r6::detail::table_7)
    {
        INFO("D= ", e.D, ", T= ", e.T, ", delta= ", e.D / critical_density, ", tau= ", critical_temperature / e.T);
        CHECK(pressure_dt (e.D, e.T) == Approx { e.P }.epsilon (1e-8));
        CHECK(massic_isochoric_heat_capacity_dt (e.D, e.T) == Approx { e.Cv }.epsilon (1e-8));
        CHECK(speed_of_sound_dt (e.D, e.T) == Approx { e.W }.epsilon (1e-8));
        CHECK(massic_entropy_dt (e.D, e.T) == Approx { e.S }.epsilon (1e-8));
    }
} // SUBCASE("main API")
SUBCASE("mixed arguments")
{
    CHECK(pressure_dt (1e3, 300.0l));
}
SUBCASE("alpha_v, kappa_T, alpha_p, beta_p")
{
    for (auto iT = 0; iT < ssize (r7::detail::T); ++iT)
    {
            const auto
        T = r7::detail::T[iT] + 273.15;
        for (auto iP = 0; iP < ssize (r7::detail::P); ++iP)
        {
                const auto
            P = r7::detail::P[iP] * 1e5;
                const auto
            D = r6_inverse::density_pt (P, T);
                const auto
            alpha_v = r7::detail::table_9[iT][iP] * 1e-6;
                const auto
            kappa_t = r7::detail::table_10[iT][iP] * 1e-9;
                const auto
            alpha_p = r7::detail::table_19[iT][iP] * 1e-3;
                const auto
            beta_p = r7::detail::table_20[iT][iP];
            INFO(
                    "T= ", T
                , ", P= ", P
                , ", D= ", D
                , ", alpha_v = ", alpha_v
                , ", kappa_t = ", kappa_t
                , ", alpha_p = ", alpha_p
                , ", beta_p = ", beta_p
            );
            CHECK(isobaric_cubic_expansion_coefficient_dt (D, T) == Approx { alpha_v }.scale (fabs (alpha_v)).epsilon (1e-2));
            CHECK(isothermal_compressibility_dt (D, T) == Approx { kappa_t }.scale (fabs (kappa_t)).epsilon (1e-2));
            CHECK(relative_pressure_coefficient_dt (D, T) == Approx { alpha_p }.scale (fabs (alpha_p)).epsilon (1e-2));
            CHECK(isothermal_stress_coefficient_dt (D, T) == Approx { beta_p }.scale (fabs (beta_p)).epsilon (1e-2));
            // Check that some basic relations hold
            CHECK(alpha_p * kappa_t / alpha_v == Approx { 1 / P });
            CHECK(beta_p * alpha_v / alpha_p == Approx { D }.scale (fabs (D)).epsilon (1e-2));
        }
    }
}
#if 0
SUBCASE("expansion, compressibility, etc.")
{
        using namespace calculisto::iapws::r7::detail;
    for (auto it = 1u; it != T.size (); ++it)
    {
            auto const
        t = T.at (it) + 273.15;
        for (auto ip = 0u; ip != P.size (); ++ip)
        {
                auto const
            p = P.at (ip) * 1e5;
                auto const
            A = table_9.at (it).at (ip) * 1e-6;
                auto const
            B = table_10.at (it).at (ip) * 1e-9;
                auto const
            C = table_19.at (it).at (ip) * 1e-3;
                auto const
            D = table_20.at (it).at (ip);
                auto const
            AA = r7::isobaric_cubic_expansion_coefficient_pt (p, t);
                auto const
            BB = r7::isothermal_compressibility_pt           (p, t);
                auto const
            CC = r7::relative_pressure_coefficient_pt        (p, t);
                auto const
            DD = r7::isothermal_stress_coefficient_pt        (p, t);
                auto const
            d7 = r7::density_pt (p, t);
                auto const
            [ d, i ] = r6_inverse::density_pt (p, t, 1e-4, info::convergence);
            if (i.converged)
            {
                {
                    INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7, ", rho = ", d, ", r7: ", AA);
                    //CHECK(isobaric_cubic_expansion_coefficient_dt (d, t) == Approx { A }.scale (fabs (A)).epsilon (1e-2));
                }{
                    INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7, ", rho = ", d, ", r7: ", BB);
                    //CHECK(isothermal_compressibility_dt           (d, t) == Approx { B }.scale (B).epsilon (1e-2));
                }{
                    INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7, ", rho = ", d, ", r7: ", CC);
                    CHECK(relative_pressure_coefficient_dt        (d, t) == Approx { C }.scale (fabs (C)).epsilon (1e-2));
                }{
                    INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7, ", rho = ", d, ", r7: ", DD);
                    CHECK(isothermal_stress_coefficient_dt        (d, t) == Approx { D }.scale (fabs (D)).epsilon (1e-2));
                }
            }
            else
            {
                MESSAGE("Failed with ""P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7);
                /*
                for (auto&& [ v, f, df ]: i.convergence)
                {
                    MESSAGE("  ", v, ", ", f, ", ", df);
                }
                */
            }
        }
    }
}
#endif
} // TEST_CASE("r6.hpp")
