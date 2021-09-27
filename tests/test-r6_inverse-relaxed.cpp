#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws::r6_inverse;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace isto::iapws::r6::r6_95_2016::detail;

TEST_CASE("r6_inverse.hpp (relaxed)")
{
    for(const auto& e: table_7)
    {
        CHECK (density_pt (e.P, e.T) == Approx { e.D }.epsilon (1e-8));
        CHECK (density_tp (e.T, e.P) == Approx { e.D }.epsilon (1e-8));
        /*
        CHECK (temperature_dp (e.D, e.P) == Approx { e.T }.epsilon (1e-8));
        CHECK (temperature_pd (e.P, e.D) == Approx { e.T }.epsilon (1e-8));
        */
    }
} // TEST_CASE("r6_inverse.hpp (relaxed)")
