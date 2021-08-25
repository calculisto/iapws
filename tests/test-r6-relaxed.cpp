#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6.hpp"
    using namespace isto::iapws::r6;
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
        CHECK(pressure (e.D, e.T) == Approx { e.P }.epsilon (1e-8));
        CHECK(massic_isochoric_heat_capacity (e.D, e.T) == Approx { e.Cv }.epsilon (1e-8));
        CHECK(speed_of_sound (e.D, e.T) == Approx { e.W }.epsilon (1e-8));
        CHECK(massic_entropy (e.D, e.T) == Approx { e.S }.epsilon (1e-8));
    }
} // SUBCASE("main API")
} // TEST_CASE("r6.hpp (relaxed)")
