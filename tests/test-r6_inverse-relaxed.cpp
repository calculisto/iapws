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
        using namespace isto::iapws;
    for(const auto& e: table_7)
    {
        // the convergence espilon is 1e-6.
        CHECK (density_pt (e.P, e.T) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        // With initial guess
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        // And convergence criterion
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T), 1e-6) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T), 1e-6) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        /*
        CHECK (temperature_dp (e.D, e.P) == Approx { e.T }.epsilon (1e-8));
        CHECK (temperature_pd (e.P, e.D) == Approx { e.T }.epsilon (1e-8));
        */
    }
} // TEST_CASE("r6_inverse.hpp (relaxed)")
