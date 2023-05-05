#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r10.hpp"
    using namespace isto::iapws::r10;

TEST_CASE("r10.hpp")
{
SUBCASE("base functions")
{
        using namespace detail;
    CHECK(g    (611.657, 273.16)     == Approx {  0.611784135000     }.scale (1)    .epsilon (1e-8));
    CHECK(g_p  (611.657, 273.16)     == Approx {  0.109085812737e-2  }.scale (1e-2) .epsilon (1e-8));
    CHECK(g_t  (611.657, 273.16)     == Approx {  0.122069433940e4   }.scale (1e4)  .epsilon (1e-8));
    CHECK(g_pp (611.657, 273.16)     == Approx { -0.128495941571e-12 }.scale (1e-12).epsilon (1e-8));
    CHECK(g_tp (611.657, 273.16)     == Approx {  0.174387964700e-6  }.scale (1e-6) .epsilon (1e-8));
    CHECK(g_tt (611.657, 273.16)     == Approx { -0.767602985875e1   }.scale (1e1)  .epsilon (1e-8));
    CHECK(g    (101325., 273.152519) == Approx {  0.101342740690e3   }.scale (1e3)  .epsilon (1e-8));
    CHECK(g_p  (101325., 273.152519) == Approx {  0.109084388214e-2  }.scale (1e-2) .epsilon (1e-8));
    CHECK(g_t  (101325., 273.152519) == Approx {  0.122076932550e4   }.scale (1e4)  .epsilon (1e-8));
    CHECK(g_pp (101325., 273.152519) == Approx { -0.128485364928e-12 }.scale (1e-12).epsilon (1e-8));
    CHECK(g_tp (101325., 273.152519) == Approx {  0.174362219972e-6  }.scale (1e-6) .epsilon (1e-8));
    CHECK(g_tt (101325., 273.152519) == Approx { -0.767598233365e1   }.scale (1e1)  .epsilon (1e-8));
    CHECK(g    (100e6, 100.)         == Approx { -0.222296513088e6   }.scale (1e6)  .epsilon (1e-8));
    CHECK(g_p  (100e6, 100.)         == Approx {  0.106193389260e-2  }.scale (1e-2) .epsilon (1e-8));
    CHECK(g_t  (100e6, 100.)         == Approx {  0.261195122589e4   }.scale (1e4)  .epsilon (1e-8));
    CHECK(g_pp (100e6, 100.)         == Approx { -0.941807981761e-13 }.scale (1e-13).epsilon (1e-8));
    CHECK(g_tp (100e6, 100.)         == Approx {  0.274505162488e-7  }.scale (1e-7) .epsilon (1e-8));
    CHECK(g_tt (100e6, 100.)         == Approx { -0.866333195517e1   }.scale (1e1)  .epsilon (1e-8));
} // SUBCASE("base functions")
SUBCASE("main API")
{
    CHECK(massic_enthalpy_pt               (611.657, 273.16    ) == Approx { -0.333444253966e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_pt       (611.657, 273.16    ) == Approx { -0.554468750000e-1  } .scale (1e-1) .epsilon (1e-8));
    CHECK(massic_internal_energy_pt        (611.657, 273.16    ) == Approx { -0.333444921197e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_pt                (611.657, 273.16    ) == Approx { -0.122069433940e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_pt (611.657, 273.16    ) == Approx {  0.209678431622e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(density_pt                       (611.657, 273.16    ) == Approx {  0.916709492200e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_pt   (611.657, 273.16    ) == Approx {  0.159863102566e-3  } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_pt          (611.657, 273.16    ) == Approx {  0.135714764659e7   } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_pt    (611.657, 273.16    ) == Approx {  0.117793449348e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_pt    (611.657, 273.16    ) == Approx {  0.114161597779e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(massic_enthalpy_pt               (101325., 273.152519) == Approx { -0.333354873637e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_pt       (101325., 273.152519) == Approx { -0.918701567000e1   } .scale (1e1)  .epsilon (1e-8));
    CHECK(massic_internal_energy_pt        (101325., 273.152519) == Approx { -0.333465403393e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_pt                (101325., 273.152519) == Approx { -0.122076932550e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_pt (101325., 273.152519) == Approx {  0.209671391024e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(density_pt                       (101325., 273.152519) == Approx {  0.916721463419e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_pt   (101325., 273.152519) == Approx {  0.159841589458e-3  } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_pt          (101325., 273.152519) == Approx {  0.135705899321e7   } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_pt    (101325., 273.152519) == Approx {  0.117785291765e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_pt    (101325., 273.152519) == Approx {  0.114154442556e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(massic_enthalpy_pt               (100e6, 100.        ) == Approx { -0.483491635676e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_pt       (100e6, 100.        ) == Approx { -0.328489902347e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_internal_energy_pt        (100e6, 100.        ) == Approx { -0.589685024936e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_pt                (100e6, 100.        ) == Approx { -0.261195122589e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_pt (100e6, 100.        ) == Approx {  0.866333195517e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(density_pt                       (100e6, 100.        ) == Approx {  0.941678203297e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_pt   (100e6, 100.        ) == Approx {  0.258495528207e-4  } .scale (1e-4) .epsilon (1e-8));
    CHECK(pressure_coefficient_pt          (100e6, 100.        ) == Approx {  0.291466166994e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_pt    (100e6, 100.        ) == Approx {  0.886880048115e-10 } .scale (1e-10).epsilon (1e-8));
    CHECK(isentropic_compressibility_pt    (100e6, 100.        ) == Approx {  0.886060982687e-10 } .scale (1e-10).epsilon (1e-8));

    CHECK(massic_enthalpy_tp               (273.16,     611.657) == Approx { -0.333444253966e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_tp       (273.16,     611.657) == Approx { -0.554468750000e-1  } .scale (1e-1) .epsilon (1e-8));
    CHECK(massic_internal_energy_tp        (273.16,     611.657) == Approx { -0.333444921197e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_tp                (273.16,     611.657) == Approx { -0.122069433940e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_tp (273.16,     611.657) == Approx {  0.209678431622e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(density_tp                       (273.16,     611.657) == Approx {  0.916709492200e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_tp   (273.16,     611.657) == Approx {  0.159863102566e-3  } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_tp          (273.16,     611.657) == Approx {  0.135714764659e7   } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_tp    (273.16,     611.657) == Approx {  0.117793449348e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_tp    (273.16,     611.657) == Approx {  0.114161597779e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(massic_enthalpy_tp               (273.152519, 101325.) == Approx { -0.333354873637e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_tp       (273.152519, 101325.) == Approx { -0.918701567000e1   } .scale (1e1)  .epsilon (1e-8));
    CHECK(massic_internal_energy_tp        (273.152519, 101325.) == Approx { -0.333465403393e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_tp                (273.152519, 101325.) == Approx { -0.122076932550e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_tp (273.152519, 101325.) == Approx {  0.209671391024e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(density_tp                       (273.152519, 101325.) == Approx {  0.916721463419e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_tp   (273.152519, 101325.) == Approx {  0.159841589458e-3  } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_tp          (273.152519, 101325.) == Approx {  0.135705899321e7   } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_tp    (273.152519, 101325.) == Approx {  0.117785291765e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_tp    (273.152519, 101325.) == Approx {  0.114154442556e-9  } .scale (1e-9) .epsilon (1e-8));
    CHECK(massic_enthalpy_tp               (100.,         100e6) == Approx { -0.483491635676e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_tp       (100.,         100e6) == Approx { -0.328489902347e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_internal_energy_tp        (100.,         100e6) == Approx { -0.589685024936e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_tp                (100.,         100e6) == Approx { -0.261195122589e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_tp (100.,         100e6) == Approx {  0.866333195517e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(density_tp                       (100.,         100e6) == Approx {  0.941678203297e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_tp   (100.,         100e6) == Approx {  0.258495528207e-4  } .scale (1e-4) .epsilon (1e-8));
    CHECK(pressure_coefficient_tp          (100.,         100e6) == Approx {  0.291466166994e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_tp    (100.,         100e6) == Approx {  0.886880048115e-10 } .scale (1e-10).epsilon (1e-8));
    CHECK(isentropic_compressibility_tp    (100.,         100e6) == Approx {  0.886060982687e-10 } .scale (1e-10).epsilon (1e-8));
} // SUBCASE("main API")
} // TEST_CASE("r10.hpp")
