#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws::r6_inverse;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace isto::iapws::r6::r6_95_2016::detail;

TEST_CASE("r6_inverse.hpp (constrained)")
{
    for(const auto& e: table_7)
    {
        // the convergence espilon is 1e-6, 1e3 is the scale.
        CHECK_U (density (pressure_t { e.P }, temperature_t { e.T }), density_t { e.D }, 1e-6 * 1e3);
        CHECK_U (density (temperature_t { e.T }, pressure_t { e.P }), density_t { e.D }, 1e-6 * 1e3);
        /*
        CHECK_U (temperature (density_t { e.D }, pressure_t { e.P }), temperature_t { e.D }, 1e-8);
        CHECK_U (temperature (pressure_t { e.P }, density_t { e.D }), temperature_t { e.D }, 1e-8);
        */
    }
} // TEST_CASE("r6_inverse.hpp (constrained)")
