#pragma once
#include <cmath>
    namespace calculisto::template_pow
{
        template <auto EXPONENT>
        constexpr auto
    pow (double value)
    {
            using std::pow;
        return pow (value, EXPONENT);
    }
} // namespace calculisto::template_pow
