#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r12.hpp"
    using namespace isto::iapws::r12;

TEST_CASE("r12.hpp (relaxed)")
{
SUBCASE("main API")
{
    CHECK_US(viscosity (temperature_t {  298.15 }, density_t {  998. }), viscosity_t {  889.735100e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (temperature_t {  298.15 }, density_t { 1200. }), viscosity_t { 1437.649467e-6 }, 1e-2, 1e-6);
    CHECK_US(viscosity (temperature_t {  373.15 }, density_t { 1000. }), viscosity_t {  307.883622e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (temperature_t {  433.15 }, density_t {    1. }), viscosity_t {   14.538324e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t {  433.15 }, density_t { 1000. }), viscosity_t {  217.685358e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (temperature_t {  873.15 }, density_t {    1. }), viscosity_t {   32.619287e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t {  873.15 }, density_t {  100. }), viscosity_t {   35.802262e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t {  873.15 }, density_t {  600. }), viscosity_t {   77.430195e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t { 1173.15 }, density_t {    1. }), viscosity_t {   44.217245e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t { 1173.15 }, density_t {  100. }), viscosity_t {   47.640433e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (temperature_t { 1173.15 }, density_t {  400. }), viscosity_t {   64.154608e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {  998. }, temperature_t {  298.15 }), viscosity_t {  889.735100e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (density_t { 1200. }, temperature_t {  298.15 }), viscosity_t { 1437.649467e-6 }, 1e-2, 1e-6);
    CHECK_US(viscosity (density_t { 1000. }, temperature_t {  373.15 }), viscosity_t {  307.883622e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (density_t {    1. }, temperature_t {  433.15 }), viscosity_t {   14.538324e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t { 1000. }, temperature_t {  433.15 }), viscosity_t {  217.685358e-6 }, 1e-3, 1e-6);
    CHECK_US(viscosity (density_t {    1. }, temperature_t {  873.15 }), viscosity_t {   32.619287e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {  100. }, temperature_t {  873.15 }), viscosity_t {   35.802262e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {  600. }, temperature_t {  873.15 }), viscosity_t {   77.430195e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {    1. }, temperature_t { 1173.15 }), viscosity_t {   44.217245e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {  100. }, temperature_t { 1173.15 }), viscosity_t {   47.640433e-6 }, 1e-4, 1e-6);
    CHECK_US(viscosity (density_t {  400. }, temperature_t { 1173.15 }), viscosity_t {   64.154608e-6 }, 1e-4, 1e-6);
    /* with near critical point correction
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  122. }), viscosity_t {   25.520677e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  222. }), viscosity_t {   31.337589e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  272. }), viscosity_t {   36.228143e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  322. }), viscosity_t {   42.961579e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  372. }), viscosity_t {   45.688204e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (temperature_t {  647.35 }, density_t {  422. }), viscosity_t {   49.436256e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  122. }, temperature_t {  647.35 }), viscosity_t {   25.520677e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  222. }, temperature_t {  647.35 }), viscosity_t {   31.337589e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  272. }, temperature_t {  647.35 }), viscosity_t {   36.228143e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  322. }, temperature_t {  647.35 }), viscosity_t {   42.961579e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  372. }, temperature_t {  647.35 }), viscosity_t {   45.688204e-6 }, 1e-5, 1e-8);
    CHECK_US(viscosity (density_t {  422. }, temperature_t {  647.35 }), viscosity_t {   49.436256e-6 }, 1e-5, 1e-8);
    */
} // SUBCASE("main API")
} // TEST_CASE("r12.hpp (relaxed)")
