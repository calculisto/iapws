#pragma once
#include "detail/common.hpp"
#include "r7.hpp"

    namespace
calculisto::iapws::lines::saturation
{
    namespace
r7 
{
    auto
pressure_t (auto const& temperature)
{
    return iapws::r7::saturation_pressure_t (temperature);
}

    constexpr auto
temperature_p (auto const& pressure)
{
    return iapws::r7::saturation_temperature_p (pressure);
}
} // namespace r7
    namespace
r6
{
    auto
pressure_t (auto const& temperature)
{
    // FIXME:
    throw not_yet_implemented_e {};
}

    constexpr auto
temperature_p (auto const& pressure)
{
    // FIXME:
    throw not_yet_implemented_e {};
}
} // namespace r6
} // namespace calculisto::iapws::lines::saturation
