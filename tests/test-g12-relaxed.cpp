#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/g12.hpp"
    using namespace isto::iapws::g12;

TEST_CASE("g12.hpp (relaxed)")
{
        using namespace detail;
    CHECK(L     (273.15 / T_LL - 1. , 0.101325e6 / rho_0 / R / T_LL) == Approx { 0.62120474 } );
    CHECK(L     (235.15 / T_LL - 1. , 0.101325e6 / rho_0 / R / T_LL) == Approx { 0.09176368 } );
    CHECK(L     (250.   / T_LL - 1.   , 200e6    / rho_0 / R / T_LL) == Approx { 0.72377081 } );
    CHECK(L     (200.   / T_LL - 1.   , 400e6    / rho_0 / R / T_LL) == Approx { 1.1553965 } );
    CHECK(L     (250.   / T_LL - 1.   , 400e6    / rho_0 / R / T_LL) == Approx { 1.4345145 } );
    CHECK(get_x (273.15 / T_LL - 1. , 0.101325e6 / rho_0 / R / T_LL) == Approx { 0.0966547155 } );
    CHECK(get_x (235.15 / T_LL - 1. , 0.101325e6 / rho_0 / R / T_LL) == Approx { 0.2551028587 } );
    CHECK(get_x (250.   / T_LL - 1.   , 200e6    / rho_0 / R / T_LL) == Approx { 0.0304292667 } );
    CHECK(get_x (200.   / T_LL - 1.   , 400e6    / rho_0 / R / T_LL) == Approx { 0.0071700809 } );
    CHECK(get_x (250.   / T_LL - 1.   , 400e6    / rho_0 / R / T_LL) == Approx { 0.0053588366 } );
    CHECK(density_tp                       ( 273.15 , 0.101325e6) == Approx {  999.84229 }.scale (1e4).epsilon (1e-6));
    CHECK(density_tp                       ( 235.15 , 0.101325e6) == Approx {  968.09999 }.scale (1e4).epsilon (1e-6));
    CHECK(density_tp                       ( 250    , 200e6)      == Approx {  1090.45677 }.scale (1e4).epsilon (1e-6));
    CHECK(density_tp                       ( 200    , 400e6)      == Approx {  1185.02800 }.scale (1e4).epsilon (1e-6));
    CHECK(density_tp                       ( 250    , 400e6)      == Approx {  1151.71517 }.scale (1e4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_tp ( 273.15 , 0.101325e6) == Approx {  -0.683042e-4 }.scale (1e-5).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_tp ( 235.15 , 0.101325e6) == Approx {  -29.63381e-4 }.scale (1e-2).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_tp ( 250    , 200e6)      == Approx {  3.267768e-4  }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_tp ( 200    , 400e6)      == Approx {  6.716009e-4  }.scale (1e-4).epsilon (1e-6));
    CHECK(thermal_expansion_coefficient_tp ( 250    , 400e6)      == Approx {  4.929927e-4  }.scale (1e-4).epsilon (1e-6));
    CHECK(isothermal_compressibility_tp    ( 273.15 , 0.101325e6) == Approx {  5.088499e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_tp    ( 235.15 , 0.101325e6) == Approx {  11.580785e-10 }.scale (1e-9).epsilon (1e-6));
    CHECK(isothermal_compressibility_tp    ( 250    , 200e6)      == Approx {  3.361311e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_tp    ( 200    , 400e6)      == Approx {  2.567237e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(isothermal_compressibility_tp    ( 250    , 400e6)      == Approx {  2.277029e-10  }.scale (1e-10).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_tp ( 273.15 , 0.101325e6) == Approx {  4218.3002  }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_tp ( 235.15 , 0.101325e6) == Approx {  5997.5632  }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_tp ( 250    , 200e6)      == Approx {  3708.3902  }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_tp ( 200    , 400e6)      == Approx {  3338.5250  }.scale (1e4).epsilon (1e-6));
    CHECK(massic_isobaric_heat_capacity_tp ( 250    , 400e6)      == Approx {  3757.2144  }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_tp                ( 273.15 , 0.101325e6) == Approx {  1402.3886  }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_tp                ( 235.15 , 0.101325e6) == Approx {  1134.5855  }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_tp                ( 250    , 200e6)      == Approx {  1668.2020  }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_tp                ( 200    , 400e6)      == Approx {  1899.3294  }.scale (1e4).epsilon (1e-6));
    CHECK(speed_of_sound_tp                ( 250    , 400e6)      == Approx {  2015.8782  }.scale (1e4).epsilon (1e-6));

    MESSAGE(homogeneous_ice_nucleation_limit_temperature_p (10e6));
    MESSAGE(homogeneous_ice_nucleation_limit_temperature_p (100e6));
    MESSAGE(homogeneous_ice_nucleation_limit_temperature_p (200e6));
    MESSAGE(homogeneous_ice_nucleation_limit_temperature_p (300e6));
}
