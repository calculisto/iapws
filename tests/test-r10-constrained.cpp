#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r10.hpp"
    using namespace isto::iapws::r10;

TEST_CASE("r10.hpp (constrained)")
{
SUBCASE("main API")
{
        using namespace isto::units::unit;
    CHECK_US(massic_enthalpy               (pressure_t { 611.657 }, temperature_t { 273.16 }    ), massic_enthalpy_t               { -0.333444253966e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (pressure_t { 611.657 }, temperature_t { 273.16 }    ), massic_energy_t                 { -0.554468750000e-1  }, 1e-1 , 1e-8);
    CHECK_US(massic_internal_energy        (pressure_t { 611.657 }, temperature_t { 273.16 }    ), massic_energy_t                 { -0.333444921197e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (pressure_t { 611.657 }, temperature_t { 273.16 }    ), massic_entropy_t                { -0.122069433940e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (pressure_t { 611.657 }, temperature_t { 273.16 }    ), massic_heat_capacity_t          {  0.209678431622e4   }, 1e4  , 1e-8);
    CHECK_US(density                       (pressure_t { 611.657 }, temperature_t { 273.16 }    ), density_t                       {  0.916709492200e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (pressure_t { 611.657 }, temperature_t { 273.16 }    ), (0.159863102566e-3 / kelvin <>)                        , 1e-3 , 1e-8);
    CHECK_US(pressure_coefficient          (pressure_t { 611.657 }, temperature_t { 273.16 }    ), (0.135714764659e7 * pascal <> / kelvin <>)             , 1e7  , 1e-8);
    CHECK_US(isothermal_compressibility    (pressure_t { 611.657 }, temperature_t { 273.16 }    ), compressibility_t               {  0.117793449348e-9  }, 1e-9 , 1e-8);
    CHECK_US(isentropic_compressibility    (pressure_t { 611.657 }, temperature_t { 273.16 }    ), compressibility_t               {  0.114161597779e-9  }, 1e-9 , 1e-8);
    CHECK_US(massic_enthalpy               (pressure_t { 101325. }, temperature_t { 273.152519 }), massic_enthalpy_t               { -0.333354873637e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (pressure_t { 101325. }, temperature_t { 273.152519 }), massic_energy_t                 { -0.918701567000e1   }, 1e1  , 1e-8);
    CHECK_US(massic_internal_energy        (pressure_t { 101325. }, temperature_t { 273.152519 }), massic_energy_t                 { -0.333465403393e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (pressure_t { 101325. }, temperature_t { 273.152519 }), massic_entropy_t                { -0.122076932550e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (pressure_t { 101325. }, temperature_t { 273.152519 }), massic_heat_capacity_t          {  0.209671391024e4   }, 1e4  , 1e-8);
    CHECK_US(density                       (pressure_t { 101325. }, temperature_t { 273.152519 }), density_t                       {  0.916721463419e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (pressure_t { 101325. }, temperature_t { 273.152519 }), (0.159841589458e-3 / kelvin <>)                        , 1e-3 , 1e-8);
    CHECK_US(pressure_coefficient          (pressure_t { 101325. }, temperature_t { 273.152519 }), (0.135705899321e7  * pascal <> / kelvin <>)            , 1e7  , 1e-8);
    CHECK_US(isothermal_compressibility    (pressure_t { 101325. }, temperature_t { 273.152519 }), compressibility_t               {  0.117785291765e-9  }, 1e-9 , 1e-8);
    CHECK_US(isentropic_compressibility    (pressure_t { 101325. }, temperature_t { 273.152519 }), compressibility_t               {  0.114154442556e-9  }, 1e-9 , 1e-8);
    CHECK_US(massic_enthalpy               (pressure_t { 100e6 }, temperature_t { 100. }        ), massic_enthalpy_t               { -0.483491635676e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (pressure_t { 100e6 }, temperature_t { 100. }        ), massic_energy_t                 { -0.328489902347e6   }, 1e6  , 1e-8);
    CHECK_US(massic_internal_energy        (pressure_t { 100e6 }, temperature_t { 100. }        ), massic_energy_t                 { -0.589685024936e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (pressure_t { 100e6 }, temperature_t { 100. }        ), massic_entropy_t                { -0.261195122589e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (pressure_t { 100e6 }, temperature_t { 100. }        ), massic_heat_capacity_t          {  0.866333195517e3   }, 1e3  , 1e-8);
    CHECK_US(density                       (pressure_t { 100e6 }, temperature_t { 100. }        ), density_t                       {  0.941678203297e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (pressure_t { 100e6 }, temperature_t { 100. }        ), (0.258495528207e-4 / kelvin <>)                        , 1e-4 , 1e-8);
    CHECK_US(pressure_coefficient          (pressure_t { 100e6 }, temperature_t { 100. }        ), (0.291466166994e6 * pascal <> / kelvin <>)             , 1e6  , 1e-8);
    CHECK_US(isothermal_compressibility    (pressure_t { 100e6 }, temperature_t { 100. }        ), compressibility_t               {  0.886880048115e-10 }, 1e-10, 1e-8);
    CHECK_US(isentropic_compressibility    (pressure_t { 100e6 }, temperature_t { 100. }        ), compressibility_t               {  0.886060982687e-10 }, 1e-10, 1e-8);
        
    CHECK_US(massic_enthalpy               (temperature_t { 273.16 },     pressure_t { 611.657 }), massic_enthalpy_t               { -0.333444253966e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (temperature_t { 273.16 },     pressure_t { 611.657 }), massic_energy_t                 { -0.554468750000e-1  }, 1e-1 , 1e-8);
    CHECK_US(massic_internal_energy        (temperature_t { 273.16 },     pressure_t { 611.657 }), massic_energy_t                 { -0.333444921197e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (temperature_t { 273.16 },     pressure_t { 611.657 }), massic_entropy_t                { -0.122069433940e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (temperature_t { 273.16 },     pressure_t { 611.657 }), massic_heat_capacity_t          {  0.209678431622e4   }, 1e4  , 1e-8);
    CHECK_US(density                       (temperature_t { 273.16 },     pressure_t { 611.657 }), density_t                       {  0.916709492200e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (temperature_t { 273.16 },     pressure_t { 611.657 }), (0.159863102566e-3 / kelvin <>)                        , 1e-3 , 1e-8);
    CHECK_US(pressure_coefficient          (temperature_t { 273.16 },     pressure_t { 611.657 }), (0.135714764659e7 * pascal <> / kelvin <>)             , 1e7  , 1e-8);
    CHECK_US(isothermal_compressibility    (temperature_t { 273.16 },     pressure_t { 611.657 }), compressibility_t               {  0.117793449348e-9  }, 1e-9 , 1e-8);
    CHECK_US(isentropic_compressibility    (temperature_t { 273.16 },     pressure_t { 611.657 }), compressibility_t               {  0.114161597779e-9  }, 1e-9 , 1e-8);
    CHECK_US(massic_enthalpy               (temperature_t { 273.152519 }, pressure_t { 101325. }), massic_enthalpy_t               { -0.333354873637e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (temperature_t { 273.152519 }, pressure_t { 101325. }), massic_energy_t                 { -0.918701567000e1   }, 1e1  , 1e-8);
    CHECK_US(massic_internal_energy        (temperature_t { 273.152519 }, pressure_t { 101325. }), massic_energy_t                 { -0.333465403393e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (temperature_t { 273.152519 }, pressure_t { 101325. }), massic_entropy_t                { -0.122076932550e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (temperature_t { 273.152519 }, pressure_t { 101325. }), massic_heat_capacity_t          {  0.209671391024e4   }, 1e4  , 1e-8);
    CHECK_US(density                       (temperature_t { 273.152519 }, pressure_t { 101325. }), density_t                       {  0.916721463419e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (temperature_t { 273.152519 }, pressure_t { 101325. }), (0.159841589458e-3 / kelvin <>)                        , 1e-3 , 1e-8);
    CHECK_US(pressure_coefficient          (temperature_t { 273.152519 }, pressure_t { 101325. }), (0.135705899321e7  * pascal <> / kelvin <>)            , 1e7  , 1e-8);
    CHECK_US(isothermal_compressibility    (temperature_t { 273.152519 }, pressure_t { 101325. }), compressibility_t               {  0.117785291765e-9  }, 1e-9 , 1e-8);
    CHECK_US(isentropic_compressibility    (temperature_t { 273.152519 }, pressure_t { 101325. }), compressibility_t               {  0.114154442556e-9  }, 1e-9 , 1e-8);
    CHECK_US(massic_enthalpy               (temperature_t { 100. },         pressure_t { 100e6 }), massic_enthalpy_t               { -0.483491635676e6   }, 1e6  , 1e-8);
    CHECK_US(massic_helmholtz_energy       (temperature_t { 100. },         pressure_t { 100e6 }), massic_energy_t                 { -0.328489902347e6   }, 1e6  , 1e-8);
    CHECK_US(massic_internal_energy        (temperature_t { 100. },         pressure_t { 100e6 }), massic_energy_t                 { -0.589685024936e6   }, 1e6  , 1e-8);
    CHECK_US(massic_entropy                (temperature_t { 100. },         pressure_t { 100e6 }), massic_entropy_t                { -0.261195122589e4   }, 1e4  , 1e-8);
    CHECK_US(massic_isobaric_heat_capacity (temperature_t { 100. },         pressure_t { 100e6 }), massic_heat_capacity_t          {  0.866333195517e3   }, 1e3  , 1e-8);
    CHECK_US(density                       (temperature_t { 100. },         pressure_t { 100e6 }), density_t                       {  0.941678203297e3   }, 1e3  , 1e-8);
    CHECK_US(cubic_expansion_coefficient   (temperature_t { 100. },         pressure_t { 100e6 }), (0.258495528207e-4 / kelvin <>)                        , 1e-4 , 1e-8);
    CHECK_US(pressure_coefficient          (temperature_t { 100. },         pressure_t { 100e6 }), (0.291466166994e6 * pascal <> / kelvin <>)             , 1e6  , 1e-8);
    CHECK_US(isothermal_compressibility    (temperature_t { 100. },         pressure_t { 100e6 }), compressibility_t               {  0.886880048115e-10 }, 1e-10, 1e-8);
    CHECK_US(isentropic_compressibility    (temperature_t { 100. },         pressure_t { 100e6 }), compressibility_t               {  0.886060982687e-10 }, 1e-10, 1e-8);
} // SUBCASE("main API")
} // TEST_CASE("r10.hpp (constrained)")
