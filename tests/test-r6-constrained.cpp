#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r6.hpp"
    using namespace isto::iapws::r6;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace r6_95_2016::detail;

TEST_CASE("r6.hpp (constrained)")
{
SUBCASE("main API")
{
    for(const auto& e: table_7)
    {
        CHECK_U(
              pressure (density_t { e.D }, temperature_t { e.T })
            , pressure_t { e.P }
            , 1e-8
        );
        CHECK_U(
              massic_isochoric_heat_capacity (
                  density_t { e.D }
                , temperature_t { e.T }
              )
            , massic_heat_capacity_t { e.Cv }
            , 1e-8
        );
        CHECK_U(
              speed_of_sound (density_t { e.D }, temperature_t { e.T })
            , velocity_t { e.W }
            , 1e-8
        );
        CHECK_U(
              massic_entropy (density_t { e.D }, temperature_t { e.T })
            , massic_entropy_t { e.S }
            , 1e-8
        );

        CHECK_U(
              pressure (temperature_t { e.T }, density_t { e.D })
            , pressure_t { e.P }
            , 1e-8
        );
        CHECK_U(
              massic_isochoric_heat_capacity (
                  temperature_t { e.T }
                , density_t { e.D }
              )
            , massic_heat_capacity_t { e.Cv }
            , 1e-8
        );
        CHECK_U(
              speed_of_sound (temperature_t { e.T }, density_t { e.D })
            , velocity_t { e.W }
            , 1e-8
        );
        CHECK_U(
              massic_entropy (temperature_t { e.T }, density_t { e.D })
            , massic_entropy_t { e.S }
            , 1e-8
        );
    }
} // SUBCASE("main API")
SUBCASE("mixed arguments")
{
    CHECK_U(pressure (density_t { 0.9965560e3 }, temperature_t { 300.0l }), pressure_t { 0.992418352e-1 * 1e6l }, 1e-8);
}
} // TEST_CASE("r6.hpp (constrained)")
