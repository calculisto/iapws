#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r8.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::iapws::r8;
    using namespace calculisto::iapws::r6_inverse;

TEST_CASE("r8.hpp")
{
SUBCASE("details")
{
        using calculisto::iapws::r8::molar_mass_of_water;
    CHECK(density_pt (   0.101325e6, 240.) == Approx { 54.33701e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt (   0.101325e6, 300.) == Approx { 55.31735e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt (  10e6       , 300.) == Approx { 55.56148e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt (1000e6       , 300.) == Approx { 68.69265e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt (  10e6       , 650.) == Approx {  2.24692e3 * molar_mass_of_water }.scale (1e1).epsilon (1e-7));
    CHECK(density_pt ( 100e6       , 650.) == Approx { 40.31090e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt ( 500e6       , 650.) == Approx { 52.58636e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK(density_pt (  10e6       , 870.) == Approx {  1.45275e3 * molar_mass_of_water }.scale (1e1).epsilon (1e-5));
    CHECK(density_pt ( 100e6       , 870.) == Approx { 20.98927e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-6));
    CHECK(density_pt ( 500e6       , 870.) == Approx { 45.01376e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
}
SUBCASE("main API")
{
    CHECK(relative_permittivity_dt (density_pt (   0.101325e6, 240.), 240.) == Approx { 104.34982 }.scale (1e2).epsilon (1e-7));
    CHECK(relative_permittivity_dt (density_pt (   0.101325e6, 300.), 300.) == Approx {  77.74735 }.scale (1e1).epsilon (1e-7));
    CHECK(relative_permittivity_dt (density_pt (  10e6       , 300.), 300.) == Approx {  78.11269 }.scale (1e1).epsilon (1e-7));
    CHECK(relative_permittivity_dt (density_pt (1000e6       , 300.), 300.) == Approx { 103.69632 }.scale (1e2).epsilon (1e-7));
    CHECK(relative_permittivity_dt (density_pt (  10e6       , 650.), 650.) == Approx {   1.26715 }.scale (1e0).epsilon (1e-5));
    CHECK(relative_permittivity_dt (density_pt ( 100e6       , 650.), 650.) == Approx {  17.71733 }.scale (1e1).epsilon (1e-6));
    CHECK(relative_permittivity_dt (density_pt ( 500e6       , 650.), 650.) == Approx {  26.62132 }.scale (1e1).epsilon (1e-7));
    CHECK(relative_permittivity_dt (density_pt (  10e6       , 870.), 870.) == Approx {   1.12721 }.scale (1e0).epsilon (1e-6));
    CHECK(relative_permittivity_dt (density_pt ( 100e6       , 870.), 870.) == Approx {   4.98281 }.scale (1e0).epsilon (1e-6));
    CHECK(relative_permittivity_dt (density_pt ( 500e6       , 870.), 870.) == Approx {  15.09746 }.scale (1e1).epsilon (1e-7));
} // SUBCASE("main API")
} // TEST_CASE("r10.hpp")


