#include <doctest/doctest.h>
#define ISTO_IAPWS_FORCE_RELAXED 1
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

TEST_CASE("r6.hpp (relaxed)")
{
SUBCASE("base functions")
{
    eq_up_to <7> (phi_0    (838.025 / 322., 647.096 / 500.),  0.204797734e1);
    eq_up_to <8> (phi_r    (838.025 / 322., 647.096 / 500.), -0.342693206e1);
    eq_up_to <8> (phi_0_d  (838.025 / 322., 647.096 / 500.),  0.384236747  );
    eq_up_to <8> (phi_r_d  (838.025 / 322., 647.096 / 500.), -0.364366650  );
    eq_up_to <9> (phi_0_dd (838.025 / 322., 647.096 / 500.), -0.147637878  );
    eq_up_to <9> (phi_r_dd (838.025 / 322., 647.096 / 500.),  0.856063701  );
    eq_up_to <8> (phi_0_t  (838.025 / 322., 647.096 / 500.),  0.904611106e1);
    eq_up_to <8> (phi_r_t  (838.025 / 322., 647.096 / 500.), -0.581403435e1);
    eq_up_to <9> (phi_0_tt (838.025 / 322., 647.096 / 500.), -0.193249185e1);
    eq_up_to <8> (phi_r_tt (838.025 / 322., 647.096 / 500.), -0.223440737e1);
    eq_up_to <9> (phi_0_dt (838.025 / 322., 647.096 / 500.),  0.           );
    eq_up_to <8> (phi_r_dt (838.025 / 322., 647.096 / 500.), -0.112176915e1);

    eq_up_to <8> (phi_0    (358.000 / 322., 647.096 / 647.), -0.156319605e1);
    eq_up_to <8> (phi_r    (358.000 / 322., 647.096 / 647.), -0.121202657e1);
    eq_up_to <9> (phi_0_d  (358.000 / 322., 647.096 / 647.),  0.899441341  );
    eq_up_to <9> (phi_r_d  (358.000 / 322., 647.096 / 647.), -0.714012024  );
    eq_up_to <9> (phi_0_dd (358.000 / 322., 647.096 / 647.), -0.808994726  );
    eq_up_to <9> (phi_r_dd (358.000 / 322., 647.096 / 647.),  0.475730696  );
    eq_up_to <8> (phi_0_t  (358.000 / 322., 647.096 / 647.),  0.980343918e1);
    eq_up_to <8> (phi_r_t  (358.000 / 322., 647.096 / 647.), -0.321722501e1);
    eq_up_to <8> (phi_0_tt (358.000 / 322., 647.096 / 647.), -0.343316334e1);
    eq_up_to <8> (phi_r_tt (358.000 / 322., 647.096 / 647.), -0.996029507e1);
    eq_up_to <9> (phi_0_dt (358.000 / 322., 647.096 / 647.),  0.           );
    eq_up_to <8> (phi_r_dt (358.000 / 322., 647.096 / 647.), -0.133214720e1);

} // SUBCASE("base functions")
SUBCASE("main API")
{
    for(const auto& e: table_7)
    {
        eq_up_to <6> (pressure (e.D, e.T), e.P);
        eq_up_to <6> (massic_isochoric_heat_capacity (e.D, e.T), e.Cv);
        eq_up_to <6> (speed_of_sound (e.D, e.T), e.W);
        eq_up_to <6> (massic_entropy (e.D, e.T), e.S);
    }
} // SUBCASE("main API")
} // TEST_CASE("r6.hpp (relaxed)")
