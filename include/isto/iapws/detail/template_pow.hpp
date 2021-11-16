#pragma once
#include <cmath>
    namespace isto::template_pow
{
        template <auto EXPONENT>
        constexpr auto
    pow (double value)
    {
            using std::pow;
        return pow (value, EXPONENT);
    }
} // namespace isto::template_pow
