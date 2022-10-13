#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/g12.hpp"
    using namespace isto::iapws::g12;

TEST_CASE("g12.hpp (constrained)")
{
    CHECK_US(density                       ( temperature_t { 273.15 } , pressure_t { 0.101325e6 }) , density_t { 999.84229 }, 1e4, 1e-6);
    CHECK_US(density                       ( temperature_t { 235.15 } , pressure_t { 0.101325e6 }) , density_t { 968.09999 }, 1e4, 1e-6);
    CHECK_US(density                       ( temperature_t { 250 }    , pressure_t { 200e6 })      , density_t { 1090.45677 }, 1e4, 1e-6);
    CHECK_US(density                       ( temperature_t { 200 }    , pressure_t { 400e6 })      , density_t { 1185.02800 }, 1e4, 1e-6);
    CHECK_US(density                       ( temperature_t { 250 }    , pressure_t { 400e6 })      , density_t { 1151.71517 }, 1e4, 1e-6);
    CHECK_US(thermal_expansion_coefficient ( temperature_t { 273.15 } , pressure_t { 0.101325e6 }) , thermal_expansion_t { -0.683042e-4 }, 1e-5, 1e-6);
    CHECK_US(thermal_expansion_coefficient ( temperature_t { 235.15 } , pressure_t { 0.101325e6 }) , thermal_expansion_t { -29.63381e-4 }, 1e-2 , 1e-6);
    CHECK_US(thermal_expansion_coefficient ( temperature_t { 250 }    , pressure_t { 200e6 })      , thermal_expansion_t { 3.267768e-4  }, 1e-4 , 1e-6);
    CHECK_US(thermal_expansion_coefficient ( temperature_t { 200 }    , pressure_t { 400e6 })      , thermal_expansion_t { 6.716009e-4  }, 1e-4 , 1e-6);
    CHECK_US(thermal_expansion_coefficient ( temperature_t { 250 }    , pressure_t { 400e6 })      , thermal_expansion_t { 4.929927e-4  }, 1e-4 , 1e-6);
    CHECK_US(isothermal_compressibility    ( temperature_t { 273.15 } , pressure_t { 0.101325e6 }) , compressibility_t { 5.088499e-10  }, 1e-10, 1e-6);
    CHECK_US(isothermal_compressibility    ( temperature_t { 235.15 } , pressure_t { 0.101325e6 }) , compressibility_t { 11.580785e-10 }, 1e-9, 1e-6);
    CHECK_US(isothermal_compressibility    ( temperature_t { 250 }    , pressure_t { 200e6 })      , compressibility_t { 3.361311e-10  }, 1e-10, 1e-6);
    CHECK_US(isothermal_compressibility    ( temperature_t { 200 }    , pressure_t { 400e6 })      , compressibility_t { 2.567237e-10  }, 1e-10, 1e-6);
    CHECK_US(isothermal_compressibility    ( temperature_t { 250 }    , pressure_t { 400e6 })      , compressibility_t { 2.277029e-10  }, 1e-10, 1e-6);
    CHECK_US(massic_isobaric_heat_capacity ( temperature_t { 273.15 } , pressure_t { 0.101325e6 }) , massic_heat_capacity_t { 4218.3002  }, 1e4, 1e-6);
    CHECK_US(massic_isobaric_heat_capacity ( temperature_t { 235.15 } , pressure_t { 0.101325e6 }) , massic_heat_capacity_t { 5997.5632  }, 1e4, 1e-6);
    CHECK_US(massic_isobaric_heat_capacity ( temperature_t { 250 }    , pressure_t { 200e6 })      , massic_heat_capacity_t { 3708.3902  }, 1e4, 1e-6);
    CHECK_US(massic_isobaric_heat_capacity ( temperature_t { 200 }    , pressure_t { 400e6 })      , massic_heat_capacity_t { 3338.5250  }, 1e4, 1e-6);
    CHECK_US(massic_isobaric_heat_capacity ( temperature_t { 250 }    , pressure_t { 400e6 })      , massic_heat_capacity_t { 3757.2144  }, 1e4, 1e-6);
    CHECK_US(speed_of_sound                ( temperature_t { 273.15 } , pressure_t { 0.101325e6 }) , velocity_t { 1402.3886  }, 1e4, 1e-6);
    CHECK_US(speed_of_sound                ( temperature_t { 235.15 } , pressure_t { 0.101325e6 }) , velocity_t { 1134.5855  }, 1e4, 1e-6);
    CHECK_US(speed_of_sound                ( temperature_t { 250 }    , pressure_t { 200e6 })      , velocity_t { 1668.2020  }, 1e4, 1e-6);
    CHECK_US(speed_of_sound                ( temperature_t { 200 }    , pressure_t { 400e6 })      , velocity_t { 1899.3294  }, 1e4, 1e-6);
    CHECK_US(speed_of_sound                ( temperature_t { 250 }    , pressure_t { 400e6 })      , velocity_t { 2015.8782  }, 1e4, 1e-6);
}
