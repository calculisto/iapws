#include <doctest/doctest.h>
#include "../include/isto/iapws/r6.hpp"
    using namespace isto::iapws::r6;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace r6_95_2016::detail;

    template <int Decimal>
    auto
eq_up_to (double a, double b)
{
        static const auto
    f = pow (1., Decimal);
    CHECK(round (f * a) == round (f * b));
}

    template <int Decimal, class T, class U>
    auto
eq_up_to_u (T const& a, U const& b)
{
    static_assert (T::dimension == U::dimension);
    eq_up_to <Decimal> (a.magnitude, b.magnitude);
}

TEST_CASE("r6.hpp (constrained)")
{
SUBCASE("main API")
{
    for(const auto& e: table_7)
    {
        eq_up_to_u <6> (
              pressure (density_t { e.D }, temperature_t { e.T })
            , pressure_t { e.P }
        );
        eq_up_to_u <6> (
              massic_isochoric_heat_capacity (
                  density_t { e.D }
                , temperature_t { e.T }
              )
            , massic_heat_capacity_t { e.Cv }
        );
        eq_up_to_u <6> (
              speed_of_sound (density_t { e.D }, temperature_t { e.T })
            , velocity_t { e.W }
        );
        eq_up_to_u <6> (
              massic_entropy (density_t { e.D }, temperature_t { e.T })
            , massic_entropy_t { e.S }
        );

        eq_up_to_u <6> (
              pressure (temperature_t { e.T }, density_t { e.D })
            , pressure_t { e.P }
        );
        eq_up_to_u <6> (
              massic_isochoric_heat_capacity (
                  temperature_t { e.T }
                , density_t { e.D }
              )
            , massic_heat_capacity_t { e.Cv }
        );
        eq_up_to_u <6> (
              speed_of_sound (temperature_t { e.T }, density_t { e.D })
            , velocity_t { e.W }
        );
        eq_up_to_u <6> (
              massic_entropy (temperature_t { e.T }, density_t { e.D })
            , massic_entropy_t { e.S }
        );
    }
} // SUBCASE("main API")
} // TEST_CASE("r6.hpp (constrained)")
