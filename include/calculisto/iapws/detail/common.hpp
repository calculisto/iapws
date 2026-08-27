#pragma once
#include "calculisto/array/array.hpp"
    using calculisto::array::array_t;

    namespace
calculisto::thermodynamics::iapws
{
    struct not_implemented_e {};
    struct not_yet_implemented_e {};
    struct internal_error_e {};

        constexpr auto
    massic_gas_constant = 0.46151805e3;

        constexpr auto
    critical_temperature = 647.096;

        constexpr auto
    critical_pressure = 22.064e6;

        constexpr auto
    critical_density = 322.0;

        constexpr auto
    triple_point_temperature = 273.16;

        constexpr auto
    triple_point_pressure = 611.657;

        constexpr auto
    boltzmann_constant = 1.380658e-23;

        constexpr auto
    avogadro_number = 6.0221367e23;

        constexpr auto
    molar_mass = 0.018015268;
    
        constexpr auto
    critical_molar_density = 322 / molar_mass;
}
