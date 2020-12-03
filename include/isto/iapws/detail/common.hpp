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

