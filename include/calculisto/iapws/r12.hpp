#pragma once
#include "r6.hpp"

#include <fmt/format.h>

    namespace 
calculisto::iapws::r12
{
    inline namespace 
r12_08_2008
{
    namespace
detail
{
        constexpr auto
    mu_0 (auto const& t)
    {
            using std::sqrt, std::pow;
        return 100. * sqrt (t) / (
              1.67752 
            + 2.20462   / t 
            + 0.6366564 / t / t 
            - 0.241605  / t / t / t
        );
    }
        constexpr auto
    mu_1 (auto t, auto d)
    {
        constexpr auto H00 =  5.20094e-1;
        constexpr auto H01 =  2.22531e-1;
        constexpr auto H02 = -2.81378e-1;
        constexpr auto H03 =  1.61913e-1;
        constexpr auto H04 = -3.25372e-2;
        
        constexpr auto H10 =  8.50895e-2;
        constexpr auto H11 =  9.99115e-1;
        constexpr auto H12 = -9.06851e-1;
        constexpr auto H13 =  2.57399e-1;
        
        constexpr auto H20 = -1.08374;
        constexpr auto H21 =  1.88797;
        constexpr auto H22 = -7.72479e-1;
        
        constexpr auto H30 = -2.89555e-1;
        constexpr auto H31 =  1.26613;
        constexpr auto H32 = -4.89837e-1;
        constexpr auto H34 =  6.98452e-2;
        constexpr auto H36 = -4.35673e-3;
        
        constexpr auto H42 = -2.57040e-1;
        constexpr auto H45 =  8.72102e-3;
        
        constexpr auto H51 =  1.20573e-1;
        constexpr auto H56 = -5.93264e-4;

            using std::exp, std::pow;
        t = (1. / t) - 1.;
        d = d - 1.;
        return exp ((d + 1.) * (
                           (H00 + H01 * d + H02 * pow (d, 2) + H03 * pow (d, 3) + H04 * pow (d, 4)                                      )
            + t          * (H10 + H11 * d + H12 * pow (d, 2) + H13 * pow (d, 3)                                                         )
            + pow (t, 2) * (H20 + H21 * d + H22 * pow (d, 2)                                                                            )
            + pow (t, 3) * (H30 + H31 * d + H32 * pow (d, 2)                    + H34 * pow (d, 4)                    + H36 * pow (d, 6))
            + pow (t, 4) * (                H42 * pow (d, 2)                                       + H45 * pow (d, 5)                   )
            + pow (t, 5) * (      H51 * d                                                                             + H56 * pow (d, 6))
        ));
    }
    
        constexpr auto
    xi (auto const& temperature, auto const& density)
    {
            constexpr auto
        nu = 0.630;
            constexpr auto
        gamma = 1.239;
            constexpr auto
        xi_0 = 0.13e-9;
            constexpr auto
        gamma_0 = 0.06;
            constexpr auto

        T_b_r = 1.5;
            auto const
        T_b = temperature / (647.096);
            auto const
        rho_b = density / (322.);
            auto const
        T_r = T_b_r * (647.096);
            auto const
        p = r6::pressure_td (temperature, density);
            auto const
        p_b = p / (22.064e6);
            auto const
        p_R = r6::pressure_td (T_r, density);
            auto const
        p_b_R = p_R / (22.064e6);
            auto const
        A = rho_b / p_b * density / r6::isothermal_stress_coefficient_dt (density, temperature);
            auto const
        B = rho_b / p_b_R * density / r6::isothermal_stress_coefficient_dt (density, T_r);
            auto const
        DX = rho_b * (A - T_b_r / T_b * B);
        return DX < 0. ? 0. : xi_0 * pow (DX / gamma_0, nu / gamma);
    }
    
        constexpr auto
    mu_2 (auto const& temperature, auto const& density)
    {
        // This requires a reduced conpressibility.
        // https://doi.org/10.1063/1.470718
        /* "A fundamental equation for the thermodynamic properties of the 
         *  fluid is needed to calculate the reduced compressibility xT* 
         *  in Eq. (6.10)"
         */
            using std::pow, std::fabs, std::tan, std::sqrt, std::log, std::atan, std::sin, std::exp;
            constexpr auto
        x_mu = 0.068;
            constexpr auto
        q_c = 1. / 1.9e-9;
            constexpr auto
        q_d = 1. / 1.1e-9;

            auto const
        x = xi (temperature, density);
            auto const
        psi_d = acos (pow (1. + q_d * q_d * x * x, -0.5));
            auto const
        w = sqrt (fabs ((q_c * x - 1.) / (q_c * x + 1.))) * tan (psi_d / 2.);
            auto const
        Lw = q_c * x > 1. ? log ((1. + w) / (1. - w)) : 2 * atan (fabs (w));
            auto const
        a = q_c * x;
            auto const
        b = q_d * x;
            auto const
        Y = x > 0.3817016416e-9 
            ? 1. / 12. * sin (3. * psi_d) - 1. / 4. / a * sin (2. * psi_d) + 1. / a / a * (1. - 5. / 4. * a * a) * sin (psi_d) - 1 / a / a / a * ((1. - 3. / 2. * a * a) * psi_d - pow (fabs (a * a - 1.), 3. / 2.) * Lw)
            : 1. / 5. * a * pow (b, 5.) * (1. - a + a * a - 765. / 504. * b * b)
        ;
        return exp (x_mu * Y);
    }
} // namespace detail
    constexpr auto
viscosity_td (auto const& temperature, auto const& density)
{
        using namespace detail;
        auto const
    t = temperature / (647.096);
        auto const
    d = density / (322.);
    return (mu_0 (t) * mu_1 (t, d) * mu_2 (temperature, density)) * (1e-6);
}
    constexpr auto
viscosity_dt (auto const& density, auto const& temperature)
{
    return viscosity_td (temperature, density);
}
} // namespace r12_08_2008
} // namespace calculisto::iapws::r12
