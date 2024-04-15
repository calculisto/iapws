#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r10_inverse.hpp"
    using namespace calculisto::iapws::r10_inverse;

TEST_CASE("r10_inverse.hpp")
{
SUBCASE("main API")
{
    CHECK(pressure_dt                      (0.916709492200e3, 273.16) == Approx {  611.657           } .scale (1e2)  .epsilon (1e-5)); // ?
    CHECK(massic_enthalpy_dt               (0.916709492200e3, 273.16) == Approx { -0.333444253966e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_dt       (0.916709492200e3, 273.16) == Approx { -0.554468750000e-1 } .scale (1e-1) .epsilon (1e-8));
    CHECK(massic_internal_energy_dt        (0.916709492200e3, 273.16) == Approx { -0.333444921197e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_dt                (0.916709492200e3, 273.16) == Approx { -0.122069433940e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_dt (0.916709492200e3, 273.16) == Approx {  0.209678431622e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_dt   (0.916709492200e3, 273.16) == Approx {  0.159863102566e-3 } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_dt          (0.916709492200e3, 273.16) == Approx {  0.135714764659e7  } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_dt    (0.916709492200e3, 273.16) == Approx {  0.117793449348e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_dt    (0.916709492200e3, 273.16) == Approx {  0.114161597779e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(pressure_dt                      (0.916721463419e3, 273.152519) == Approx {  101325.           } .scale (1e5)  .epsilon (1e-8));
    CHECK(massic_enthalpy_dt               (0.916721463419e3, 273.152519) == Approx { -0.333354873637e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_dt       (0.916721463419e3, 273.152519) == Approx { -0.918701567000e1  } .scale (1e1)  .epsilon (1e-8));
    CHECK(massic_internal_energy_dt        (0.916721463419e3, 273.152519) == Approx { -0.333465403393e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_dt                (0.916721463419e3, 273.152519) == Approx { -0.122076932550e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_dt (0.916721463419e3, 273.152519) == Approx {  0.209671391024e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_dt   (0.916721463419e3, 273.152519) == Approx {  0.159841589458e-3 } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_dt          (0.916721463419e3, 273.152519) == Approx {  0.135705899321e7  } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_dt    (0.916721463419e3, 273.152519) == Approx {  0.117785291765e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_dt    (0.916721463419e3, 273.152519) == Approx {  0.114154442556e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(pressure_dt                      (0.941678203297e3, 100.) == Approx {  100e6              } .scale (1e8)  .epsilon (1e-8));
    CHECK(massic_enthalpy_dt               (0.941678203297e3, 100.) == Approx { -0.483491635676e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_dt       (0.941678203297e3, 100.) == Approx { -0.328489902347e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_internal_energy_dt        (0.941678203297e3, 100.) == Approx { -0.589685024936e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_dt                (0.941678203297e3, 100.) == Approx { -0.261195122589e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_dt (0.941678203297e3, 100.) == Approx {  0.866333195517e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_dt   (0.941678203297e3, 100.) == Approx {  0.258495528207e-4  } .scale (1e-4) .epsilon (1e-8));
    CHECK(pressure_coefficient_dt          (0.941678203297e3, 100.) == Approx {  0.291466166994e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_dt    (0.941678203297e3, 100.) == Approx {  0.886880048115e-10 } .scale (1e-10).epsilon (1e-8));
    CHECK(isentropic_compressibility_dt    (0.941678203297e3, 100.) == Approx {  0.886060982687e-10 } .scale (1e-10).epsilon (1e-8));

    CHECK(pressure_td                      (273.16, 0.916709492200e3) == Approx {  611.657           } .scale (1e2)  .epsilon (1e-5)); // ?
    CHECK(massic_enthalpy_td               (273.16, 0.916709492200e3) == Approx { -0.333444253966e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_td       (273.16, 0.916709492200e3) == Approx { -0.554468750000e-1 } .scale (1e-1) .epsilon (1e-8));
    CHECK(massic_internal_energy_td        (273.16, 0.916709492200e3) == Approx { -0.333444921197e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_td                (273.16, 0.916709492200e3) == Approx { -0.122069433940e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_td (273.16, 0.916709492200e3) == Approx {  0.209678431622e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_td   (273.16, 0.916709492200e3) == Approx {  0.159863102566e-3 } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_td          (273.16, 0.916709492200e3) == Approx {  0.135714764659e7  } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_td    (273.16, 0.916709492200e3) == Approx {  0.117793449348e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_td    (273.16, 0.916709492200e3) == Approx {  0.114161597779e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(pressure_td                      (273.152519, 0.916721463419e3) == Approx {  101325.           } .scale (1e5)  .epsilon (1e-8));
    CHECK(massic_enthalpy_td               (273.152519, 0.916721463419e3) == Approx { -0.333354873637e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_td       (273.152519, 0.916721463419e3) == Approx { -0.918701567000e1  } .scale (1e1)  .epsilon (1e-8));
    CHECK(massic_internal_energy_td        (273.152519, 0.916721463419e3) == Approx { -0.333465403393e6  } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_td                (273.152519, 0.916721463419e3) == Approx { -0.122076932550e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_td (273.152519, 0.916721463419e3) == Approx {  0.209671391024e4  } .scale (1e4)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_td   (273.152519, 0.916721463419e3) == Approx {  0.159841589458e-3 } .scale (1e-3) .epsilon (1e-8));
    CHECK(pressure_coefficient_td          (273.152519, 0.916721463419e3) == Approx {  0.135705899321e7  } .scale (1e7)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_td    (273.152519, 0.916721463419e3) == Approx {  0.117785291765e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(isentropic_compressibility_td    (273.152519, 0.916721463419e3) == Approx {  0.114154442556e-9 } .scale (1e-9) .epsilon (1e-8));
    CHECK(pressure_td                      (100., 0.941678203297e3) == Approx {  100e6              } .scale (1e8)  .epsilon (1e-8));
    CHECK(massic_enthalpy_td               (100., 0.941678203297e3) == Approx { -0.483491635676e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_helmholtz_energy_td       (100., 0.941678203297e3) == Approx { -0.328489902347e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_internal_energy_td        (100., 0.941678203297e3) == Approx { -0.589685024936e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(massic_entropy_td                (100., 0.941678203297e3) == Approx { -0.261195122589e4   } .scale (1e4)  .epsilon (1e-8));
    CHECK(massic_isobaric_heat_capacity_td (100., 0.941678203297e3) == Approx {  0.866333195517e3   } .scale (1e3)  .epsilon (1e-8));
    CHECK(cubic_expansion_coefficient_td   (100., 0.941678203297e3) == Approx {  0.258495528207e-4  } .scale (1e-4) .epsilon (1e-8));
    CHECK(pressure_coefficient_td          (100., 0.941678203297e3) == Approx {  0.291466166994e6   } .scale (1e6)  .epsilon (1e-8));
    CHECK(isothermal_compressibility_td    (100., 0.941678203297e3) == Approx {  0.886880048115e-10 } .scale (1e-10).epsilon (1e-8));
    CHECK(isentropic_compressibility_td    (100., 0.941678203297e3) == Approx {  0.886060982687e-10 } .scale (1e-10).epsilon (1e-8));
} // SUBCASE("main API")
} // TEST_CASE("r10.hpp")
