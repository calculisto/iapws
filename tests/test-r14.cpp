#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r14.hpp"
    using namespace calculisto::iapws::r14;

TEST_CASE("r14.hpp")
{
SUBCASE("main API")
{
    CHECK( ih::melting_pressure_t (260.0) == Approx {  138.268e6 }.scale (1e8).epsilon (1e-6));
    CHECK(iii::melting_pressure_t (254.0) == Approx {  268.685e6 }.scale (1e8).epsilon (1e-6));
    CHECK(  v::melting_pressure_t (265.0) == Approx {  479.640e6 }.scale (1e8).epsilon (1e-6));
    CHECK( vi::melting_pressure_t (320.0) == Approx { 1356.760e6 }.scale (1e9).epsilon (1e-5));
    CHECK(vii::melting_pressure_t (550.0) == Approx { 6308.710e6 }.scale (1e9).epsilon (1e-6));
    CHECK( ih::sublimation_pressure_t (230.0) == Approx { 8.94735 }.scale (1e0).epsilon (1e-6));

    CHECK(iii::melting_temperature_p ( 268.685e6) == Approx { 254.0 }.scale (1e3).epsilon (1e-6));
    CHECK(  v::melting_temperature_p ( 479.640e6) == Approx { 265.0 }.scale (1e8).epsilon (1e-6));
    CHECK( vi::melting_temperature_p (1356.760e6) == Approx { 320.0 }.scale (1e9).epsilon (1e-5));

/* Not possible.
    CHECK(melting_pressure_t (260.0) == Approx {  138.268e6 }.scale (1e8).epsilon (1e-8));
    CHECK(melting_pressure_t (254.0) == Approx {  268.685e6 }.scale (1e8).epsilon (1e-8));
    CHECK(melting_pressure_t (265.0) == Approx {  479.640e6 }.scale (1e8).epsilon (1e-8));
    CHECK(melting_pressure_t (320.0) == Approx { 1356.760e6 }.scale (1e9).epsilon (1e-8));
    CHECK(melting_pressure_t (550.0) == Approx { 6305.710e6 }.scale (1e9).epsilon (1e-8));
    CHECK(sublimation_pressure_t (230.0) == Approx { 8.94735e-6 }.scale (1e-6).epsilon (1e-8));
*/
} // SUBCASE("main API")
} // TEST_CASE("r14.hpp")
