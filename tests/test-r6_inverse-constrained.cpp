#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws::r6_inverse;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace isto::iapws::r6::r6_95_2016::detail;

TEST_CASE("r6_inverse.hpp (constrained)")
{
        using namespace isto::iapws;
    for(const auto& e: table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T);
        // the convergence espilon is 1e-6, 1e3 is the scale.
        CHECK_U (density (pressure_t { e.P }, temperature_t { e.T }), density_t { e.D }, 1e-6 * 1e3);
        CHECK_U (density (temperature_t { e.T }, pressure_t { e.P }), density_t { e.D }, 1e-6 * 1e3);
        // with initial_guess
        CHECK_U (density (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), density_t { e.D }, 1e3 * 1e-6);
        CHECK_U (density (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), density_t { e.D }, 1e3 * 1e-6);
        // with convergence criterion
        CHECK_U (density (pressure_t { e.P }, temperature_t { e.T }), density_t { e.D }, 1e3 * 1e-6);
        CHECK_U (density (temperature_t { e.T }, pressure_t { e.P }), density_t { e.D }, 1e3 * 1e-6);
        // with both
        CHECK_U (density (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), density_t { e.D }, 1e3 * 1e-6);
        CHECK_U (density (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), density_t { e.D }, 1e3 * 1e-6);
    }
    {
            const auto
        [ r, i ] = density_tp (temperature_t { 300.0 },  pressure_t { 0.992418352e-1 * 1e6 }, info::convergence);
        CHECK(i.convergence.size () > 1);
        /*
        for (auto&& [ v, f, df ]: i.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
        */
    }
    for(const auto& e: table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T);
        CHECK_U (massic_isochoric_heat_capacity (pressure_t { e.P }, temperature_t { e.T }), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (temperature_t { e.T }, pressure_t { e.P }), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);
        CHECK_U (massic_isochoric_heat_capacity (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, massic_heat_capacity_t { e.Cv }, 1e3 * 1e-6);

        CHECK_U (speed_of_sound (pressure_t { e.P }, temperature_t { e.T }), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (temperature_t { e.T }, pressure_t { e.P }), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, velocity_t { e.W }, 1e3 * 1e-6);
        CHECK_U (speed_of_sound (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, velocity_t { e.W }, 1e3 * 1e-6);

        CHECK_U (massic_entropy (pressure_t { e.P }, temperature_t { e.T }), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (temperature_t { e.T }, pressure_t { e.P }), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T })), massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (pressure_t { e.P }, temperature_t { e.T }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, massic_entropy_t { e.S }, 1e3 * 1e-6);
        CHECK_U (massic_entropy (temperature_t { e.T }, pressure_t { e.P }, r7::density (pressure_t { e.P }, temperature_t { e.T }), info::iterations).first, massic_entropy_t { e.S }, 1e3 * 1e-6);

    }
} // TEST_CASE("r6_inverse.hpp (constrained)")
