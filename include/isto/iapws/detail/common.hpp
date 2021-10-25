#pragma once
#include <isto/array/array.hpp>
    using isto::array::array_t;
#if __has_include(<isto/units/units.hpp>) && !defined ISTO_IAPWS_FORCE_RELAXED
#   include <isto/units/units.hpp>
        using namespace isto::units;
#   define ISTO_IAPWS_FLAVOR_CONSTRAINED 1
#   undef ISTO_IAPWS_FLAVOR_RELAXED 
#else
#   define ISTO_IAPWS_FLAVOR_RELAXED 1
#   undef ISTO_IAPWS_FLAVOR_CONSTRAINED 
#endif

    struct not_implemented_e {};
    struct not_yet_implemented_e {};
    struct internal_error_e {};

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_U_GC * unit::joule <> / unit::kelvin <> / unit::kilogram <>
#define ISTO_IAPWS_U_T  * unit::kelvin <>
#define ISTO_IAPWS_U_P  * unit::pascal <>
#define ISTO_IAPWS_U_D  * unit::kilogram <> / pow <3> (unit::metre <>)
#define ISTO_IAPWS_U_V  * pow <3> (unit::metre <>) / unit::kilogram <>
#define ISTO_IAPWS_U_H  * unit::joule <> / unit::kilogram <>
#define ISTO_IAPWS_U_S  * unit::joule <> / unit::kilogram <> / unit::kelvin <>
#define ISTO_IAPWS_U_M  * unit::pascal <> * unit::second <>
#else
#define ISTO_IAPWS_U_GC
#define ISTO_IAPWS_U_T
#define ISTO_IAPWS_U_P
#define ISTO_IAPWS_U_D
#define ISTO_IAPWS_U_V
#define ISTO_IAPWS_U_H
#define ISTO_IAPWS_U_S
#define ISTO_IAPWS_U_M
#endif

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED

#define ISTO_IAPWS_P1 pressure_t <T>
#define ISTO_IAPWS_T1 temperature_t <T>
#define ISTO_IAPWS_D1 density_t <T>
#define ISTO_IAPWS_H1 massic_enthalpy_t <T>
#define ISTO_IAPWS_S1 massic_entropy_t <T>

#define ISTO_IAPWS_P2 pressure_t <U>
#define ISTO_IAPWS_T2 temperature_t <U>
#define ISTO_IAPWS_D2 density_t <U>
#define ISTO_IAPWS_H2 massic_enthalpy_t <U>
#define ISTO_IAPWS_S2 massic_entropy_t <U>

#define ISTO_IAPWS_P3 pressure_t <V>
#define ISTO_IAPWS_T3 temperature_t <V>
#define ISTO_IAPWS_D3 density_t <V>
#define ISTO_IAPWS_H3 massic_enthalpy_t <V>
#define ISTO_IAPWS_S3 massic_entropy_t <V>

#define ISTO_IAPWS_P4 pressure_t <W>
#define ISTO_IAPWS_T4 temperature_t <W>
#define ISTO_IAPWS_D4 density_t <W>
#define ISTO_IAPWS_H4 massic_enthalpy_t <W>
#define ISTO_IAPWS_S4 massic_entropy_t <W>

#else

#define ISTO_IAPWS_P1 T
#define ISTO_IAPWS_T1 T
#define ISTO_IAPWS_D1 T
#define ISTO_IAPWS_H1 T
#define ISTO_IAPWS_S1 T

#define ISTO_IAPWS_P2 U
#define ISTO_IAPWS_T2 U
#define ISTO_IAPWS_D2 U
#define ISTO_IAPWS_H2 U
#define ISTO_IAPWS_S2 U

#define ISTO_IAPWS_P3 V
#define ISTO_IAPWS_T3 V
#define ISTO_IAPWS_D3 V
#define ISTO_IAPWS_H3 V
#define ISTO_IAPWS_S3 V

#define ISTO_IAPWS_P4 W
#define ISTO_IAPWS_T4 W
#define ISTO_IAPWS_D4 W
#define ISTO_IAPWS_H4 W
#define ISTO_IAPWS_S4 W

#endif

