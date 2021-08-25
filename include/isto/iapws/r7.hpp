#pragma once
#include "detail/common.hpp"

    namespace 
isto::iapws::r7
{
    inline namespace 
r7_97_2012
{

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_U_GC * unit::joule <> / unit::kelvin <> / unit::kilogram <>
#define ISTO_IAPWS_U_T  * unit::kelvin <>
#define ISTO_IAPWS_U_P  * unit::pascal <>
#define ISTO_IAPWS_U_D  * unit::kilogram <> / pow <3> (unit::metre <>)
#else
#define ISTO_IAPWS_U_GC
#define ISTO_IAPWS_U_T
#define ISTO_IAPWS_U_P
#define ISTO_IAPWS_U_D
#endif

    constexpr auto
massic_gas_constant = 0.461526 ISTO_IAPWS_U_GC;

    constexpr auto
critical_temperature = 647.096 ISTO_IAPWS_U_T;

    constexpr auto
critical_pressure = 22.064e6 ISTO_IAPWS_U_P;

    constexpr auto
critical_density = 322.0 ISTO_IAPWS_U_D;

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_PRESSURE    pressure_t <T>
#define ISTO_IAPWS_TEMPERATURE temperature_t <T>
#define ISTO_IAPWS_DENSITY     density_t <T>
#else
#define ISTO_IAPWS_PRESSURE    T
#define ISTO_IAPWS_TEMPERATURE T
#define ISTO_IAPWS_DENSITY     T
#endif

    namespace
detail
{

// §4 Auxiliary Equation for the Boundary between Regions 2 and 3.
    template <class T>  
    constexpr auto
b23 (ISTO_IAPWS_TEMPERATURE const& temperature)
{
        constexpr double
    n1 =  0.34805185628969e3;
        constexpr double 
    n2 = -0.11671859879975e1;
        constexpr double 
    n3 =  0.10192970039326e-2;
        auto const
    t = temperature / (1 ISTO_IAPWS_U_T);
    return (n1 + n2 * t + n3 * t * t) * (1e6 ISTO_IAPWS_U_P);
}

    template <class T>
    /*constexpr*/ auto // Someday, std::sqrt might be constexpr.
b23i (ISTO_IAPWS_PRESSURE const& pressure)
{
        using std::sqrt;
        constexpr double 
    n3 = 0.10192970039326e-2;
        constexpr double 
    n4 = 0.57254459862746e3;
        constexpr double 
    n5 = 0.13918839778870e2;
        auto const
    p = pressure / (1e6 ISTO_IAPWS_U_P);
    return (n4 + sqrt ((p - n5) / n3)) * (1 ISTO_IAPWS_U_T);
}

    namespace
r4
{
    namespace
detail
{
    constexpr auto
n = array_t //{{{
{
       0.11670521452767e4  
    , -0.72421316703206e6  
    , -0.17073846940092e2  
    ,  0.12020824702470e5  
    , -0.32325550322333e7  
    ,  0.14915108613530e2 
    , -0.48232657361591e4 
    ,  0.40511340542057e6 
    , -0.23855557567849 
    ,  0.65017534844798e3 
}; //}}}
} // namespace detail
} // namespace r4
} // namespace detail

    template <class T>
    constexpr auto
saturation_pressure (ISTO_IAPWS_TEMPERATURE const& temperature)
{
        using std::pow;
        using std::sqrt;
        using namespace detail::r4::detail;
        auto const
    t = temperature / (1 ISTO_IAPWS_U_T);
        auto const
    theta = t + n[8] / (t - n[9]);
        auto const
    A = theta * theta + n[0] * theta + n[1];
        auto const
    B = n[2] * theta * theta + n[3] * theta + n[4];
        auto const
    C = n[5] * theta * theta + n[6] * theta + n[7];
    return pow (2. * C / (-B + sqrt (B * B - 4. * A * C)), 4) * 1e6 ISTO_IAPWS_U_P;
}

    template <class T>
    constexpr auto
saturation_temperature (ISTO_IAPWS_PRESSURE const& pressure)
{
        using std::pow;
        using std::sqrt;
        using namespace detail::r4::detail;
        auto const
    beta = pow (pressure / (1e6 ISTO_IAPWS_U_P), 0.25);
        auto const
    E = beta * beta + n[2] * beta + n[5];
        auto const
    F = n[0] * beta * beta + n[3] * beta + n[6];
        auto const
    G = n[1] * beta * beta + n[4] * beta + n[7];
        auto const
    D = 2. * G / (-F - sqrt (F * F - 4. * E * G));
    return (n[9] + D - sqrt (pow (n[9] + D, 2) - 4. * (n[8] + n[9] * D))) / 2. ISTO_IAPWS_U_T;
}

    template <class T>
    constexpr int
region (ISTO_IAPWS_PRESSURE const& pressure, ISTO_IAPWS_TEMPERATURE const& temperature)
{
    if (temperature <= (623.15 ISTO_IAPWS_U_T))
    {
            auto const
        ps = saturation_pressure (temperature);
        if (pressure > ps) return 1;
        if (pressure < ps) return 2;
        return 4;
    }
    if (temperature < 1073.15 ISTO_IAPWS_U_T)
    {
            auto const
        p23 = detail::b23 (temperature);
        if (pressure < p23) return 2;
        return 3;
    }
    return 5;
}

#if 0

#define ISTO_THERMODYNAMICS_IAPWSR7_GEN(NAME)                                                           \
    template <class T>                                                                                  \
    auto                                                                                                \
NAME (ISTO_IAPWS_PRESSURE const& pressure, ISTO_IAPWS_TEMPERATURE const& temperature)                   \
{                                                                                                       \
        using namespace detail;                                                                         \
    switch (region (pressure, temperature))                                                             \
    {                                                                                                   \
        case 1: return r1::NAME (pressure, temperature);                                                \
        case 2: return r2::NAME (pressure, temperature);                                                \
        case 3: throw not_implemented_e { /*"You are in region 3, inversion needed"*/ };                \
        case 4: throw not_implemented_e { /*"On the saturation line, this is probably unreachable"*/ }; \
        case 5: return r5::NAME (pressure, temperature);                                                \
        default: throw not_yet_implemented_e {};                                                        \
    }                                                                                                   \
}                                                                                                       \
    template <class T>                                                                                  \
    auto                                                                                                \
NAME (ISTO_IAPWS_TEMPERATURE const& temperature, ISTO_IAPWS_PRESSURE const& pressure)                   \
{                                                                                                       \
    return NAME (pressure, temperature);                                                                \
}                                                                                                       \

ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_volume)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_enthalpy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_internal_energy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_entropy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_isobaric_heat_capacity)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_isochoric_heat_capacity)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(speed_of_sound)

#undef ISTO_THERMODYNAMICS_IAPWSR7_GEN

    template <class T>
    auto
density (ISTO_IAPWS_PRESSURE const& pressure, ISTO_IAPWS_TEMPERATURE const& temperature)
{
    return pow <-1> (massic_volume (pressure, temperature));
}
    template <class T>
    auto
density (ISTO_IAPWS_TEMPERATURE const& temperature, ISTO_IAPWS_PRESSURE const& pressure)
{
    return density (pressure, temperature);
}

#endif

#undef ISTO_IAPWS_PRESSURE
#undef ISTO_IAPWS_TEMPERATURE
#undef ISTO_IAPWS_DENSITY

#undef ISTO_IAPWS_U_GC
#undef ISTO_IAPWS_U_T
#undef ISTO_IAPWS_U_P
#undef ISTO_IAPWS_U_D

} // inline namespace r7_97_2012
} // namespace isto::iapws::r7
// vim: foldmethod=marker
