#pragma once
#include <numbers>
#include "detail/common.hpp"

    namespace
calculisto::iapws::r8
{
    inline namespace
r8_09_1997
{
    namespace
detail
{
        constexpr auto
    N = array_t
    {
           0.978224486826
        , -0.957771379375
        ,  0.237511794148
        ,  0.714692244396
        , -0.298217036956
        , -0.108863472196
        ,  0.949327488264e-1
        , -0.980469816509e-2
        ,  0.165167634970e-4
        ,  0.937359795772e-4
        , -0.123179218720e-9
    };

        constexpr auto
    N_12 = 0.196096504426e-2;

        constexpr auto
    i = array_t
    {
          1.
        , 1.
        , 1.
        , 2.
        , 3.
        , 3.
        , 4.
        , 5.
        , 6.
        , 7.
        , 10.
    };

        constexpr auto
    j = array_t
    {
          0.25
        , 1.
        , 2.5
        , 1.5
        , 1.5
        , 2.5
        , 2.
        , 2.
        , 5.
        , 0.5
        , 10.
    };

        constexpr auto
    permittivity_of_free_space = 
        1. / (4e-7 * std::numbers::pi * 299792458 * 299792458);

        constexpr auto
    mean_molecular_polarizability = 1.636e-40;

        constexpr auto
    molecular_dipole_moment = 6.138e-30;

        constexpr auto
    boltzmann_constant = 1.380658e-23;

        constexpr auto
    avogadro_number = 6.0221367e23;

        constexpr auto
    molar_mass_of_water = 0.018015268;
    
        constexpr auto
    critical_molar_density = 322 / molar_mass_of_water;

        constexpr auto
    critical_temperature = 647.096;

        constexpr auto
    g (auto const& molar_density, auto const& temperature)
    {
            using std::pow;
            const auto
        rho = molar_density / critical_molar_density;
            const auto
        tau = critical_temperature / temperature;
        return 
              1 
            + sum (N * pow (rho, i) * pow (tau, j)) 
            + N_12 * rho * pow (temperature / 228 - 1, -1.2)
        ;
    }
        constexpr auto
    A (auto const& molar_density, auto const& temperature)
    {
        return 
              avogadro_number 
            * molecular_dipole_moment 
            * molecular_dipole_moment 
            * molar_density 
            * g (molar_density, temperature) 
            / permittivity_of_free_space 
            / boltzmann_constant 
            / temperature
        ;
    }
        constexpr auto
    B (auto const& molar_density)
    {
        return
              avogadro_number
            * mean_molecular_polarizability
            * molar_density
            / 3
            / permittivity_of_free_space
        ;
    }
        constexpr auto
    epsilon (auto const& molar_density, auto const& temperature)
    {
            const auto
        a = A (molar_density, temperature);
            const auto
        b = B (molar_density);
            using std::sqrt;
        return 
            (1 + a + 5 * b + sqrt (
                9 + 2 * a + 18 * b + a * a + 10 * a * b + 9 * b * b
            )) / (4 - 4 * b)
        ;
    }
} // namespace detail

    using namespace detail;

    constexpr auto
relative_permittivity_dt (auto const& density, auto const& temperature)
{
    return epsilon (density / molar_mass_of_water, temperature);
}

} // namespace r8_09_1997
} // namespace calculisto::iapws::r8
