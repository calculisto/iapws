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
#   define ISTO_IAPWS_U_GC * unit::joule <> / unit::kelvin <> / unit::kilogram <>
#   define ISTO_IAPWS_U_T  * unit::kelvin <>
#   define ISTO_IAPWS_U_P  * unit::pascal <>
#   define ISTO_IAPWS_U_D  * unit::kilogram <> / pow <3> (unit::metre <>)
#   define ISTO_IAPWS_U_V  * pow <3> (unit::metre <>) / unit::kilogram <>
#   define ISTO_IAPWS_U_H  * unit::joule <> / unit::kilogram <>
#   define ISTO_IAPWS_U_S  * unit::joule <> / unit::kilogram <> / unit::kelvin <>
#   define ISTO_IAPWS_U_M  * unit::pascal <> * unit::second <>
#else
#   define ISTO_IAPWS_U_GC
#   define ISTO_IAPWS_U_T
#   define ISTO_IAPWS_U_P
#   define ISTO_IAPWS_U_D
#   define ISTO_IAPWS_U_V
#   define ISTO_IAPWS_U_H
#   define ISTO_IAPWS_U_S
#   define ISTO_IAPWS_U_M
#endif

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#   define ISTO_IAPWS_P Pressure
#   define ISTO_IAPWS_T Temperature
#   define ISTO_IAPWS_D Density
#   define ISTO_IAPWS_H Massic_enthalpy
#   define ISTO_IAPWS_S Massic_entropy
#else
#   define ISTO_IAPWS_P
#   define ISTO_IAPWS_T
#   define ISTO_IAPWS_D
#   define ISTO_IAPWS_H
#   define ISTO_IAPWS_S
#endif

