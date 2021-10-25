#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws;
    using namespace isto::iapws::r6;
#include "../include/isto/iapws/r7.hpp"
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace r6_95_2016::detail;

TEST_CASE("r6.hpp (relaxed)")
{
SUBCASE("base functions")
{
    CHECK(phi_0    (838.025 / 322., 647.096 / 500.) == Approx {  0.204797734e1 }.epsilon (1e-8));
    CHECK(phi_r    (838.025 / 322., 647.096 / 500.) == Approx { -0.342693206e1 }.epsilon (1e-8));
    CHECK(phi_0_d  (838.025 / 322., 647.096 / 500.) == Approx {  0.384236747   }.epsilon (1e-8));
    CHECK(phi_r_d  (838.025 / 322., 647.096 / 500.) == Approx { -0.364366650   }.epsilon (1e-8));
    CHECK(phi_0_dd (838.025 / 322., 647.096 / 500.) == Approx { -0.147637878   }.epsilon (1e-8));
    CHECK(phi_r_dd (838.025 / 322., 647.096 / 500.) == Approx {  0.856063701   }.epsilon (1e-8));
    CHECK(phi_0_t  (838.025 / 322., 647.096 / 500.) == Approx {  0.904611106e1 }.epsilon (1e-8));
    CHECK(phi_r_t  (838.025 / 322., 647.096 / 500.) == Approx { -0.581403435e1 }.epsilon (1e-8));
    CHECK(phi_0_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.193249185e1 }.epsilon (1e-8));
    CHECK(phi_r_tt (838.025 / 322., 647.096 / 500.) == Approx { -0.223440737e1 }.epsilon (1e-8));
    CHECK(phi_0_dt (838.025 / 322., 647.096 / 500.) == Approx {  0.            }.epsilon (1e-8));
    CHECK(phi_r_dt (838.025 / 322., 647.096 / 500.) == Approx { -0.112176915e1 }.epsilon (1e-8));
    CHECK(phi_0    (358.000 / 322., 647.096 / 647.) == Approx { -0.156319605e1 }.epsilon (1e-8));
    CHECK(phi_r    (358.000 / 322., 647.096 / 647.) == Approx { -0.121202657e1 }.epsilon (1e-8));
    CHECK(phi_0_d  (358.000 / 322., 647.096 / 647.) == Approx {  0.899441341   }.epsilon (1e-8));
    CHECK(phi_r_d  (358.000 / 322., 647.096 / 647.) == Approx { -0.714012024   }.epsilon (1e-8));
    CHECK(phi_0_dd (358.000 / 322., 647.096 / 647.) == Approx { -0.808994726   }.epsilon (1e-8));
    CHECK(phi_r_dd (358.000 / 322., 647.096 / 647.) == Approx {  0.475730696   }.epsilon (1e-8));
    CHECK(phi_0_t  (358.000 / 322., 647.096 / 647.) == Approx {  0.980343918e1 }.epsilon (1e-8));
    CHECK(phi_r_t  (358.000 / 322., 647.096 / 647.) == Approx { -0.321722501e1 }.epsilon (1e-8));
    CHECK(phi_0_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.343316334e1 }.epsilon (1e-8));
    CHECK(phi_r_tt (358.000 / 322., 647.096 / 647.) == Approx { -0.996029507e1 }.epsilon (1e-8));
    CHECK(phi_0_dt (358.000 / 322., 647.096 / 647.) == Approx {  0.            }.epsilon (1e-8));
    CHECK(phi_r_dt (358.000 / 322., 647.096 / 647.) == Approx { -0.133214720e1 }.epsilon (1e-8));

} // SUBCASE("base functions")
SUBCASE("main API")
{
    for(const auto& e: table_7)
    {
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
SUBCASE("expansion, compressibility, etc.")
{
        using namespace isto::iapws::r7::detail;
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
                    CHECK(isobaric_cubic_expansion_coefficient_dt (d, t) == Approx { A }.scale (fabs (A)).epsilon (1e-2));
                }{
                    INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", rho_0 = ", d7, ", rho = ", d, ", r7: ", BB);
                    CHECK(isothermal_compressibility_dt           (d, t) == Approx { B }.scale (B).epsilon (1e-2));
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
} // TEST_CASE("r6.hpp (relaxed)")
