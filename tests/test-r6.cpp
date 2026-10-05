#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r6.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::thermodynamics::iapws;
    using namespace calculisto::thermodynamics::iapws::r6;
#include "../include/calculisto/iapws/detail/data_for_the_tests.hpp"
#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;

TEST_CASE("r6.hpp")
{
SUBCASE("base functions")
{
            using namespace r6::detail;
    CHECK(phi_0    (838.025 / 322., 647.096 / 500.) == Approx {  0.204797734e1 });
    CHECK(phi_r    (838.025 / 322., 647.096 / 500.) == Approx { -0.342693206e1 });
    CHECK(phi_0_d  (838.025 / 322., 647.096 / 500.) == Approx {  0.384236747   });
    CHECK(phi_r_d  (838.025 / 322., 647.096 / 500.) == Approx { -0.364366650   });
    CHECK(phi_0_dd (838.025 / 322., 647.096 / 500.) == Approx { -0.147637878   });
    CHECK(phi_r_dd (838.025 / 322., 647.096 / 500.) == Approx {  0.856063701   });
    CHECK(phi_0_t  (838.025 / 322., 647.096 / 500.) == Approx {  0.904611106e1 });
    CHECK(phi_r_t  (838.025 / 322., 647.096 / 500.) == Approx { -0.581403435e1 });
    CHECK(phi_0_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.193249185e1 });
    CHECK(phi_r_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.223440737e1 });
    //CHECK(phi_0_dt (838.025 / 322., 647.096 / 500.) == Approx {  0.           8));
    CHECK(phi_r_dt (838.025 / 322., 647.096 / 500.) == Approx { -0.112176915e1 });
    CHECK(phi_0    (358.000 / 322., 647.096 / 647.) == Approx { -0.156319605e1 });
    CHECK(phi_r    (358.000 / 322., 647.096 / 647.) == Approx { -0.121202657e1 });
    CHECK(phi_0_d  (358.000 / 322., 647.096 / 647.) == Approx {  0.899441341   });
    CHECK(phi_r_d  (358.000 / 322., 647.096 / 647.) == Approx { -0.714012024   });
    CHECK(phi_0_dd (358.000 / 322., 647.096 / 647.) == Approx { -0.808994726   });
    CHECK(phi_r_dd (358.000 / 322., 647.096 / 647.) == Approx {  0.475730696   });
    CHECK(phi_0_t  (358.000 / 322., 647.096 / 647.) == Approx {  0.980343918e1 });
    CHECK(phi_r_t  (358.000 / 322., 647.096 / 647.) == Approx { -0.321722501e1 });
    CHECK(phi_0_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.343316334e1 });
    CHECK(phi_r_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.996029507e1 });
    //CHECK(phi_0_dt (358.000 / 322., 647.096 / 647.) == Approx {  0.            });
    CHECK(phi_r_dt (358.000 / 322., 647.096 / 647.) == Approx { -0.133214720e1 });

} // SUBCASE("base functions")
SUBCASE("derivatives, including the third-orders")
{
    for(const auto& e: r6::detail::table_7)
    {
            const auto
        delta = e.D / critical_density;
            const auto
        tau = critical_temperature / e.T;
        INFO("D= ", e.D, ", T= ", e.T, ", delta= ", delta, ", tau= ", tau);

            using namespace r6::detail;
        CHECK(phi_1_r_d   (delta, tau) == Approx { central_finite_difference <0> (phi_1_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_dd  (delta, tau) == Approx { central_finite_difference <0> (phi_1_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_ddd (delta, tau) == Approx { central_finite_difference <0> (phi_1_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_t   (delta, tau) == Approx { central_finite_difference <1> (phi_1_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_tt  (delta, tau) == Approx { central_finite_difference <1> (phi_1_r_t  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_ttt (delta, tau) == Approx { central_finite_difference <1> (phi_1_r_tt <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_dt  (delta, tau) == Approx { central_finite_difference <1> (phi_1_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_ddt (delta, tau) == Approx { central_finite_difference <1> (phi_1_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_1_r_dtt (delta, tau) == Approx { central_finite_difference <0> (phi_1_r_tt <double, double>, 1e-6, delta, tau) });

        CHECK(phi_2_r_d   (delta, tau) == Approx { central_finite_difference <0> (phi_2_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_dd  (delta, tau) == Approx { central_finite_difference <0> (phi_2_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_ddd (delta, tau) == Approx { central_finite_difference <0> (phi_2_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_t   (delta, tau) == Approx { central_finite_difference <1> (phi_2_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_tt  (delta, tau) == Approx { central_finite_difference <1> (phi_2_r_t  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_ttt (delta, tau) == Approx { central_finite_difference <1> (phi_2_r_tt <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_dt  (delta, tau) == Approx { central_finite_difference <1> (phi_2_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_ddt (delta, tau) == Approx { central_finite_difference <1> (phi_2_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_2_r_dtt (delta, tau) == Approx { central_finite_difference <0> (phi_2_r_tt <double, double>, 1e-6, delta, tau) });

        CHECK(phi_3_r_d   (delta, tau) == Approx { central_finite_difference <0> (phi_3_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_dd  (delta, tau) == Approx { central_finite_difference <0> (phi_3_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_ddd (delta, tau) == Approx { central_finite_difference <0> (phi_3_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_t   (delta, tau) == Approx { central_finite_difference <1> (phi_3_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_tt  (delta, tau) == Approx { central_finite_difference <1> (phi_3_r_t  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_ttt (delta, tau) == Approx { central_finite_difference <1> (phi_3_r_tt <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_dt  (delta, tau) == Approx { central_finite_difference <1> (phi_3_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_ddt (delta, tau) == Approx { central_finite_difference <1> (phi_3_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_3_r_dtt (delta, tau) == Approx { central_finite_difference <0> (phi_3_r_tt <double, double>, 1e-6, delta, tau) });

        CHECK(sum(Delta_d   (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta    <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_dd  (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_d  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_ddd (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_dd <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_t   (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta    <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_tt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_t  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_ttt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_tt <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_dt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_d  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_ddt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_dd <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Delta_dtt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_dt <double, double>, 1e-6, delta, tau)) });

        CHECK(sum(Psi_d   (delta, tau)) == Approx { sum(central_finite_difference <0> (Psi    <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_dd  (delta, tau)) == Approx { sum(central_finite_difference <0> (Psi_d  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_ddd (delta, tau)) == Approx { sum(central_finite_difference <0> (Psi_dd <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_t   (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi    <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_tt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi_t  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_ttt (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi_tt <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_dt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi_d  <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_ddt (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi_dd <double, double>, 1e-6, delta, tau)) });
        CHECK(sum(Psi_dtt (delta, tau)) == Approx { sum(central_finite_difference <1> (Psi_dt <double, double>, 1e-6, delta, tau)) });

        CHECK(sum(Delta_b_d   (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b    <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_dd  (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b_d  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_ddd (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b_dd <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_t   (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b    <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_tt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_t  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_ttt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_tt <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_dt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_d  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_ddt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_dd <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_dtt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_dt <double, double>, 1e-7, delta, tau)) });

        CHECK(sum(Delta_b_Psi_d   (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b_Psi    <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_dd  (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b_Psi_d  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_ddd (delta, tau)) == Approx { sum(central_finite_difference <0> (Delta_b_Psi_dd <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_t   (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi    <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_tt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi_t  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_ttt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi_tt <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_dt  (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi_d  <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_ddt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi_dd <double, double>, 1e-7, delta, tau)) });
        CHECK(sum(Delta_b_Psi_dtt (delta, tau)) == Approx { sum(central_finite_difference <1> (Delta_b_Psi_dt <double, double>, 1e-7, delta, tau)) });

        CHECK(phi_4_r_d   (delta, tau) == Approx { central_finite_difference <0> (phi_4_r    <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_dd  (delta, tau) == Approx { central_finite_difference <0> (phi_4_r_d  <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_ddd (delta, tau) == Approx { central_finite_difference <0> (phi_4_r_dd <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_t   (delta, tau) == Approx { central_finite_difference <1> (phi_4_r    <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_tt  (delta, tau) == Approx { central_finite_difference <1> (phi_4_r_t  <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_ttt (delta, tau) == Approx { central_finite_difference <1> (phi_4_r_tt <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_dt  (delta, tau) == Approx { central_finite_difference <1> (phi_4_r_d  <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_ddt (delta, tau) == Approx { central_finite_difference <1> (phi_4_r_dd <double, double>, 1e-7, delta, tau) });
        CHECK(phi_4_r_dtt (delta, tau) == Approx { central_finite_difference <0> (phi_4_r_tt <double, double>, 1e-7, delta, tau) });

        CHECK(phi_r_d   (delta, tau) == Approx { central_finite_difference <0> (phi_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_dd  (delta, tau) == Approx { central_finite_difference <0> (phi_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_ddd (delta, tau) == Approx { central_finite_difference <0> (phi_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_t   (delta, tau) == Approx { central_finite_difference <1> (phi_r    <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_tt  (delta, tau) == Approx { central_finite_difference <1> (phi_r_t  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_ttt (delta, tau) == Approx { central_finite_difference <1> (phi_r_tt <double, double>, 1e-7, delta, tau) });
        CHECK(phi_r_dt  (delta, tau) == Approx { central_finite_difference <1> (phi_r_d  <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_ddt (delta, tau) == Approx { central_finite_difference <1> (phi_r_dd <double, double>, 1e-6, delta, tau) });
        CHECK(phi_r_dtt (delta, tau) == Approx { central_finite_difference <0> (phi_r_tt <double, double>, 1e-6, delta, tau) });
    }

}
SUBCASE("main API")
{
    for(const auto& e: r6::detail::table_7)
    {
        INFO("D= ", e.D, ", T= ", e.T, ", delta= ", e.D / critical_density, ", tau= ", critical_temperature / e.T);
        CHECK(pressure_dt (e.D, e.T) == Approx { e.P });
        CHECK(massic_isochoric_heat_capacity_dt (e.D, e.T) == Approx { e.Cv });
        CHECK(speed_of_sound_dt (e.D, e.T) == Approx { e.W });
        CHECK(massic_entropy_dt (e.D, e.T) == Approx { e.S });
    }
} // SUBCASE("main API")
SUBCASE("mixed arguments")
{
    CHECK(pressure_dt (1e3, 300.0l));
}
/* FIXME:: do not test that against r7
SUBCASE("alpha_v, kappa_T, alpha_p, beta_p")
{
    for (auto iT = 0u; iT < ssize (r7::detail::T); ++iT)
    {
            const auto
        T = r7::detail::T[iT] + 273.15; // °C -> k
        for (auto iP = 0u; iP < ssize (r7::detail::P); ++iP)
        {
                const auto
            P = r7::detail::P[iP] * 1e5; // bar -> Pa
                const auto
            D = r6_inverse::density_pt (P, T);
                const auto
            alpha_v_r7 = r7::detail::table_9[iT][iP] * 1e-6;
                const auto
            kappa_t_r7 = r7::detail::table_10[iT][iP] * 1e-9;
                const auto
            alpha_p_r7 = r7::detail::table_19[iT][iP] * 1e-3;
                const auto
            beta_p_r7 = r7::detail::table_20[iT][iP];
                const auto
            d_D_d_T_ft = central_finite_difference <1> (
                  [](auto P, auto T)
                  {
                    return r6_inverse::density_pt (P, T);
                  }
                , 1e-6
                , P
                , T
            );
                const auto
            alpha_v_ft = -d_D_d_T_ft / D;
                const auto
            d_D_d_P_ft = central_finite_difference <0> (
                  [](auto P, auto T)
                  {
                    return r6_inverse::density_pt (P, T);
                  }
                , 1e-6
                , P
                , T
            );
                const auto
            kappa_t_ft = d_D_d_P_ft / D;

            INFO(
                    "T= ", T
                , ", P= ", P
                , ", D= ", D
                , ", alpha_v = ", alpha_v_r7
                , ", kappa_t = ", kappa_t_r7
                , ", alpha_p = ", alpha_p_r7
                , ", beta_p = ", beta_p_r7
            );

            CHECK(isobaric_cubic_expansion_coefficient_dt (D, T) == Approx { alpha_v_r7 }.scale (fabs (alpha_v_r7)).epsilon (1e-2));
            CHECK(isothermal_compressibility_dt (D, T) == Approx { kappa_t_r7 }.scale (fabs (kappa_t_r7)).epsilon (1e-2));
            CHECK(relative_pressure_coefficient_dt (D, T) == Approx { alpha_p_r7 }.scale (fabs (alpha_p_r7)).epsilon (1e-2));
            CHECK(isothermal_stress_coefficient_dt (D, T) == Approx { beta_p_r7 }.scale (fabs (beta_p_r7)).epsilon (1e-2));
            // Check that some basic relations hold
            CHECK(alpha_p_r7 * kappa_t_r7 / alpha_v_r7 == Approx { 1 / P });
            CHECK(beta_p_r7 * alpha_v_r7 / alpha_p_r7 == Approx { D }.scale (fabs (D)).epsilon (1e-2));
        }
    }
}
*/
SUBCASE("Derivatives needed by the Born functions")
{
    for(const auto& e: r6::detail::table_7)
    {
        INFO("D= ", e.D, ", T= ", e.T, ", P= ", e.P);
            const auto
        d_P_d_D_at_T = d_pressure_d_density_at_temperature_dt (e.D, e.T);
            const auto
        d_P_d_D_at_T_fd = central_finite_difference (
              r6::pressure_dt <double, double>
            , 1e-5
            , e.D
            , e.T
        );
        CHECK(d_P_d_D_at_T == Approx { d_P_d_D_at_T_fd  });

        // NOTE: d_density_d_temperature_at_pressure_dt is tested in r6_inverse.cpp
        // NOTE: d2_density_d2_temperature_at_pressure_dt is tested in r6_inverse.cpp
        
            const auto
        d_D_d_TT_at_PD = d_density_d2_temperature_at_pressure_density_dt (e.D, e.T);
            const auto
        d_D_d_TT_at_PD_fd = central_finite_difference <1> (
              d_density_d_temperature_at_pressure_dt <double, double>
            , 1e-6
            , e.D
            , e.T
        );
        CHECK(d_D_d_TT_at_PD == Approx { d_D_d_TT_at_PD_fd });
        {
                const auto
            delta = e.D / critical_density;
                const auto
            tau = critical_temperature / e.T;
            CHECK(r6::detail::d_dddt_numerator_d_d (delta, tau) == Approx { central_finite_difference <0> (
                  r6::detail::dddt_numerator <double, double>
                , 1e-6
                , delta
                , tau
            ) });
            CHECK(r6::detail::d_dddt_denominator_d_d (delta, tau) == Approx { central_finite_difference <0> (
                  r6::detail::dddt_denominator <double, double>
                , 1e-6
                , delta
                , tau
            ) });
        }

            const auto
        d_D_d_TD_at_PT = d_density_d_temperature_density_at_pressure_temperature_dt (e.D, e.T);
            const auto
        d_D_d_TD_at_PT_fd = central_finite_difference <0> (
              r6::d_density_d_temperature_at_pressure_dt <double, double>
            , 1e-6
            , e.D
            , e.T
        );
        CHECK(d_D_d_TD_at_PT == Approx { d_D_d_TD_at_PT_fd  });
    }
}
#if 0
SUBCASE("expansion, compressibility, etc.")
{
        using namespace calculisto::thermodynamics::iapws::r7::detail;
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
