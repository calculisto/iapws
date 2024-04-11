#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/g12_inverse.hpp"
    using namespace isto::iapws::g12_inverse;

TEST_CASE("g12.hpp")
{
    CHECK(pressure_td                      ( 273.15 , 999.84229 ) == Approx {  0.101325e6    }.scale (1e5).epsilon (1e-6));
    CHECK(pressure_td                      ( 235.15 , 968.09999 ) == Approx {  0.101325e6    }.scale (1e5).epsilon (1e-6));
    CHECK(pressure_td                      ( 250.   , 1090.45677) == Approx {  200e6         }.scale (1e8).epsilon (1e-6));
    CHECK(pressure_td                      ( 200.   , 1185.02800) == Approx {  400e6         }.scale (1e8).epsilon (1e-6));
    CHECK(pressure_td                      ( 250.   , 1151.71517) == Approx {  400e6         }.scale (1e8).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_td ( 273.15 , 999.84229 ) == Approx {  -0.683042e-4  }.scale (1e-5).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_td ( 235.15 , 968.09999 ) == Approx {  -29.63381e-4  }.scale (1e-2).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_td ( 250.   , 1090.45677) == Approx {  3.267768e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_td ( 200.   , 1185.02800) == Approx {  6.716009e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_td ( 250.   , 1151.71517) == Approx {  4.929927e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(isothermal_compressibility_td    ( 273.15 , 999.84229 ) == Approx {  5.088499e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_td    ( 235.15 , 968.09999 ) == Approx {  11.580785e-10 }.scale (1e-9).epsilon (1e-6));
    CHECK(isothermal_compressibility_td    ( 250.   , 1090.45677) == Approx {  3.361311e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_td    ( 200.   , 1185.02800) == Approx {  2.567237e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_td    ( 250.   , 1151.71517) == Approx {  2.277029e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_td ( 273.15 , 999.84229 ) == Approx {  4218.3002     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_td ( 235.15 , 968.09999 ) == Approx {  5997.5632     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_td ( 250.   , 1090.45677) == Approx {  3708.3902     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_td ( 200.   , 1185.02800) == Approx {  3338.5250     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_td ( 250.   , 1151.71517) == Approx {  3757.2144     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_td                ( 273.15 , 999.84229 ) == Approx {  1402.3886     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_td                ( 235.15 , 968.09999 ) == Approx {  1134.5855     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_td                ( 250.   , 1090.45677) == Approx {  1668.2020     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_td                ( 200.   , 1185.02800) == Approx {  1899.3294     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_td                ( 250.   , 1151.71517) == Approx {  2015.8782     }.scale (1e4).epsilon (1e-6));

    CHECK(pressure_dt                      ( 999.84229  , 273.15 ) == Approx {  0.101325e6    }.scale (1e5).epsilon (1e-6));
    CHECK(pressure_dt                      ( 968.09999  , 235.15 ) == Approx {  0.101325e6    }.scale (1e5).epsilon (1e-6));
    CHECK(pressure_dt                      ( 1090.45677 , 250.   ) == Approx {  200e6         }.scale (1e8).epsilon (1e-6));
    CHECK(pressure_dt                      ( 1185.02800 , 200.   ) == Approx {  400e6         }.scale (1e8).epsilon (1e-6));
    CHECK(pressure_dt                      ( 1151.71517 , 250.   ) == Approx {  400e6         }.scale (1e8).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_dt ( 999.84229  , 273.15 ) == Approx {  -0.683042e-4  }.scale (1e-5).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_dt ( 968.09999  , 235.15 ) == Approx {  -29.63381e-4  }.scale (1e-2).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_dt ( 1090.45677 , 250.   ) == Approx {  3.267768e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_dt ( 1185.02800 , 200.   ) == Approx {  6.716009e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_dt ( 1151.71517 , 250.   ) == Approx {  4.929927e-4   }.scale (1e-4).epsilon (1e-6));
    CHECK(isothermal_compressibility_dt    ( 999.84229  , 273.15 ) == Approx {  5.088499e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_dt    ( 968.09999  , 235.15 ) == Approx {  11.580785e-10 }.scale (1e-9).epsilon (1e-6));
    CHECK(isothermal_compressibility_dt    ( 1090.45677 , 250.   ) == Approx {  3.361311e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_dt    ( 1185.02800 , 200.   ) == Approx {  2.567237e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_dt    ( 1151.71517 , 250.   ) == Approx {  2.277029e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_dt ( 999.84229  , 273.15 ) == Approx {  4218.3002     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_dt ( 968.09999  , 235.15 ) == Approx {  5997.5632     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_dt ( 1090.45677 , 250.   ) == Approx {  3708.3902     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_dt ( 1185.02800 , 200.   ) == Approx {  3338.5250     }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_dt ( 1151.71517 , 250.   ) == Approx {  3757.2144     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_dt                ( 999.84229  , 273.15 ) == Approx {  1402.3886     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_dt                ( 968.09999  , 235.15 ) == Approx {  1134.5855     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_dt                ( 1090.45677 , 250.   ) == Approx {  1668.2020     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_dt                ( 1185.02800 , 200.   ) == Approx {  1899.3294     }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_dt                ( 1151.71517 , 250.   ) == Approx {  2015.8782     }.scale (1e4).epsilon (1e-6));
}
