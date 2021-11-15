#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r12.hpp"
    using namespace isto::iapws::r12;

TEST_CASE("r12.hpp (relaxed)")
{
SUBCASE("base functions")
{
        using namespace detail;
    CHECK(xi (647.35, 122.) == Approx {  0.309247e-9 }.scale (1e-9).epsilon (1e-6));
    CHECK(xi (647.35, 222.) == Approx {  1.571405e-9 }.scale (1e-9).epsilon (1e-7));
    CHECK(xi (647.35, 272.) == Approx {  5.266522e-9 }.scale (1e-9).epsilon (1e-5));
    CHECK(xi (647.35, 322.) == Approx { 16.590209e-9 }.scale (1e-8).epsilon (1e-8));
    CHECK(xi (647.35, 372.) == Approx {  5.603768e-9 }.scale (1e-9).epsilon (1e-5));
    CHECK(xi (647.35, 422.) == Approx {  1.876244e-9 }.scale (1e-9).epsilon (1e-5));
    CHECK(mu_2 (647.35, 122.) == Approx { 1.00000289 }.epsilon (1e-8));
    CHECK(mu_2 (647.35, 222.) == Approx { 1.00375120 }.epsilon (1e-8));
    CHECK(mu_2 (647.35, 272.) == Approx { 1.03416789 }.epsilon (1e-8));
    CHECK(mu_2 (647.35, 322.) == Approx { 1.09190440 }.epsilon (1e-8));
    CHECK(mu_2 (647.35, 372.) == Approx { 1.03665871 }.epsilon (1e-8));
    CHECK(mu_2 (647.35, 422.) == Approx { 1.00596332 }.epsilon (1e-8));
} // SUBCASE("base functions")
SUBCASE("main API")
{
    CHECK(viscosity_td ( 298.15,  998.) == Approx {  889.735100e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_td ( 298.15, 1200.) == Approx { 1437.649467e-6 }.scale (1e-2).epsilon (1e-6));
    CHECK(viscosity_td ( 373.15, 1000.) == Approx {  307.883622e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_td ( 433.15,    1.) == Approx {   14.538324e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td ( 433.15, 1000.) == Approx {  217.685358e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_td ( 873.15,    1.) == Approx {   32.619287e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td ( 873.15,  100.) == Approx {   35.802262e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td ( 873.15,  600.) == Approx {   77.430195e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td (1173.15,    1.) == Approx {   44.217245e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td (1173.15,  100.) == Approx {   47.640433e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_td (1173.15,  400.) == Approx {   64.154608e-6 }.scale (1e-4).epsilon (1e-6));

    CHECK(viscosity_dt ( 998.,  298.15) == Approx {  889.735100e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_dt (1200.,  298.15) == Approx { 1437.649467e-6 }.scale (1e-2).epsilon (1e-6));
    CHECK(viscosity_dt (1000.,  373.15) == Approx {  307.883622e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_dt (   1.,  433.15) == Approx {   14.538324e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt (1000.,  433.15) == Approx {  217.685358e-6 }.scale (1e-3).epsilon (1e-6));
    CHECK(viscosity_dt (   1.,  873.15) == Approx {   32.619287e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt ( 100.,  873.15) == Approx {   35.802262e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt ( 600.,  873.15) == Approx {   77.430195e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt (   1., 1173.15) == Approx {   44.217245e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt ( 100., 1173.15) == Approx {   47.640433e-6 }.scale (1e-4).epsilon (1e-6));
    CHECK(viscosity_dt ( 400., 1173.15) == Approx {   64.154608e-6 }.scale (1e-4).epsilon (1e-6));

    CHECK(viscosity_td (647.35, 122.) == Approx { 25.520677e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_td (647.35, 222.) == Approx { 31.337589e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_td (647.35, 272.) == Approx { 36.228143e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_td (647.35, 322.) == Approx { 42.961579e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_td (647.35, 372.) == Approx { 45.688204e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_td (647.35, 422.) == Approx { 49.436256e-6 }.scale (1e-5).epsilon (1e-8));

    CHECK(viscosity_dt (122., 647.35) == Approx { 25.520677e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_dt (222., 647.35) == Approx { 31.337589e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_dt (272., 647.35) == Approx { 36.228143e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_dt (322., 647.35) == Approx { 42.961579e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_dt (372., 647.35) == Approx { 45.688204e-6 }.scale (1e-5).epsilon (1e-8));
    CHECK(viscosity_dt (422., 647.35) == Approx { 49.436256e-6 }.scale (1e-5).epsilon (1e-8));
} // SUBCASE("main API")
} // TEST_CASE("r12.hpp (relaxed)")
