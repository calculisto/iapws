#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6_inverse.hpp"
#include "../include/isto/iapws/r7.hpp"
#include "../include/isto/iapws/r10.hpp"
#include "../include/isto/iapws/r14.hpp"
    using namespace isto::iapws;
#include <optional>
#include <map>
#include <thread>
#include <fstream>
#include <fmt/format.h>
    using fmt::print, fmt::format;
    using namespace std::literals;



    auto
diagram_r6_inv_ext (int n)
{
    {
            auto
        o = std::ofstream { "density_r6_inv_ext_ylin.dat" };
            auto
        e = std::ofstream { "density_r6_inv_ext_ylin_failed.dat" };
            auto
        c = std::ofstream { "density_r6_inv_ext_ylin_failed_convergence.txt" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r6_inv_ext_ylin.dat" };
            auto
        og  = std::ofstream { "massic_gibbs_free_energy_r6_inv_ext_ylin.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r6_inv_ext_ylin.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r6_inv_ext_ylin.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r6_inv_ext_ylin.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r6_inv_ext_ylin.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r6_inv_ext_ylin.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r6_inv_ext_ylin.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r6_inv_ext_ylin.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r6_inv_ext_ylin.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r6_inv_ext_ylin.dat" };
            const auto
        t0 = 173.15;
            const auto
        t1 = 5073.15;
            const auto
        dt = (t1 - t0) / (n - 1);
            const auto
        p0 = 0;
            const auto
        p1 = 100e9;
            const auto
        dp = (p1 - p0) / n;
        for (auto it = 0; it < n; ++it)
        {
                const auto
            t = t0 + it * dt;
                const auto
            pmih = r14::ih::melting_pressure_t (t);
                const auto
            pmiii = r14::iii::melting_pressure_t (t);
                const auto
            pmv = r14::v::melting_pressure_t (t);
                const auto
            pmvi = r14::vi::melting_pressure_t (t);
                [[ maybe_unused ]]const auto
            pmvii = r14::vii::melting_pressure_t (t);
                const auto
            ps = r14::sublimation_pressure_t (t);
                auto
            d = 1000.;
            for (auto ip = 1; ip < n; ++ip)
            {
                    const auto
                p = p0 + ip * dp;
                if (p <  611.657   && t <= 273.16 && p > ps) continue;
                if (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)  continue;
                if (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p > pmiii) continue;
                if (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p > pmv)   continue;
                if (p >= 632.4e6   && p <  2216e6    && t < 355.    && p > pmvi)  continue;
                if (p >= 2216e6    && p <  2.06e10   && t < 715.    && p > pmvii) continue;
                if (p >= 2.06e10   && t < 715.) continue;
                if (p > 100e6 || (t > 1073.15 && p > 50e6))
                {
                        const auto
                    [ r, i ] = r6_inverse::density_pt (p, t, d, p * 1e-7, info::convergence);
                    if (i.converged)
                    {
                        d = r;
                        o << t << " " << p << " " << d << "\n";
                        oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                        ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                        og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                        os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                        ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                        ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                        ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                        oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                        okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                        oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                        obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                        continue;
                    }
                    e << t << " " << p << " " << 1 << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                    const auto
                [ r, i ] = r6_inverse::density_pt (p, t, r7::density_pt (p, t), p * 1e-7, info::convergence);
                if (i.converged)
                {
                    d = r;
                    o << t << " " << p << " " << d << "\n";
                    oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                    ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                    og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                    os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                    ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                    ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                    ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                    oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                    okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                    oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                    obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                    continue;
                }
                e << t << " " << p << " " << 1 << "\n";
                c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                for (auto&& [ v, f ,d ]: i.convergence)
                {
                    c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                }
                continue;
            }
        }
    }{
            auto
        o = std::ofstream { "density_r6_inv_ext_ylog.dat" };
            auto
        e = std::ofstream { "density_r6_inv_ext_ylog_failed.dat" };
            auto
        c = std::ofstream { "density_r6_inv_ext_ylog_failed_convergence.txt" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r6_inv_ext_ylog.dat" };
            auto
        og  = std::ofstream { "massic_gibbs_free_energy_r6_inv_ext_ylog.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r6_inv_ext_ylog.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r6_inv_ext_ylog.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r6_inv_ext_ylog.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r6_inv_ext_ylog.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r6_inv_ext_ylog.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r6_inv_ext_ylog.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r6_inv_ext_ylog.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r6_inv_ext_ylog.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r6_inv_ext_ylog.dat" };
            const auto
        t0 = 173.15;
            const auto
        t1 = 5073.15;
            const auto
        dt = (t1 - t0) / (n - 1);
            const auto
        p0 = 1;
            const auto
        p1 = 100e9;
            const auto
        dp = exp (log (p1 / p0) / (n - 1));
        for (auto it = 0; it < n; ++it)
        {
                const auto
            t = t0 + it * dt;
                const auto
            pmih = r14::ih::melting_pressure_t (t);
                const auto
            pmiii = r14::iii::melting_pressure_t (t);
                const auto
            pmv = r14::v::melting_pressure_t (t);
                const auto
            pmvi = r14::vi::melting_pressure_t (t);
                [[ maybe_unused ]]const auto
            pmvii = r14::vii::melting_pressure_t (t);
                const auto
            ps = r14::sublimation_pressure_t (t);
                auto
            d = 1000.;
            for (auto ip = 0; ip < n; ++ip)
            {
                    const auto
                p = p0 * pow (dp, ip);
                if (p <  611.657   && t <= 273.16 && p > ps) continue;
                if (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)  continue;
                if (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p > pmiii) continue;
                if (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p > pmv)   continue;
                if (p >= 632.4e6   && p <  2216e6    && t < 355.    && p > pmvi)  continue;
                if (p >= 2216e6    && p <  2.06e10   && t < 715.    && p > pmvii) continue;
                if (p >= 2.06e10   && t < 715.) continue;
                if (p > 100e6 || (t > 1073.15 && p > 50e6))
                {
                        const auto
                    [r, i] = r6_inverse::density_pt (p, t, d, p * 1e-7, info::convergence);
                    if (i.converged)
                    {
                        d = r;
                        o << t << " " << p << " " << d << "\n";
                        oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                        ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                        og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                        os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                        ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                        ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                        ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                        oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                        okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                        oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                        obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                        continue;
                    }
                    e << t << " " << p << " " << 1 << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                    const auto
                [ r, i ] = r6_inverse::density_pt (p, t, r7::density_pt (p, t), p * 1e-7, info::convergence);
                if (i.converged)
                {
                    d = r;
                    o << t << " " << p << " " << d << "\n";
                    oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                    ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                    og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                    os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                    ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                    ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                    ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                    oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                    okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                    oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                    obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                    continue;
                }
                e << t << " " << p << " " << 1 << "\n";
                c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                for (auto&& [ v, f ,d ]: i.convergence)
                {
                    c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                }
                continue;
            }
        }
    }
}
    auto
diagram_r6_inv (int n)
{
    {
            auto
        o   = std::ofstream { "density_r6_inv_ylin.dat" };
            auto
        e   = std::ofstream { "density_r6_inv_ylin_failed.dat" };
            auto
        c   = std::ofstream { "density_r6_inv_ylin_failed_convergence.txt" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r6_inv_ylin.dat" };
            auto
        og  = std::ofstream { "massic_gibbs_free_energy_r6_inv_ylin.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r6_inv_ylin.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r6_inv_ylin.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r6_inv_ylin.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r6_inv_ylin.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r6_inv_ylin.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r6_inv_ylin.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r6_inv_ylin.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r6_inv_ylin.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r6_inv_ylin.dat" };

            const auto
        t0 = 173.15;
            const auto
        t1 = 1273.15;
            const auto
        dt = (t1 - t0) / (n - 1);
            const auto
        p0 = 0;
            const auto
        p1 = 1000e6;
            const auto
        dp = (p1 - p0) / (n - 1);
        for (auto it = 0; it < n; ++it)
        {
                const auto
            t = t0 + it * dt;
                const auto
            pmih = r14::ih::melting_pressure_t (t);
                const auto
            pmiii = r14::iii::melting_pressure_t (t);
                const auto
            pmv = r14::v::melting_pressure_t (t);
                const auto
            pmvi = r14::vi::melting_pressure_t (t);
                [[ maybe_unused ]]const auto
            pmvii = r14::vii::melting_pressure_t (t);
                const auto
            ps = r14::sublimation_pressure_t (t);
                auto
            d = 1000.;
            for (auto ip = 1; ip < n; ++ip)
            {
                    const auto
                p = p0 + ip * dp;
                if (p <  611.657   && t <= 273.16 && p > ps) continue;
                if (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)  continue;
                if (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p > pmiii) continue;
                if (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p > pmv)   continue;
                if (p >= 632.4e6   && p <  2216e6    && t < 355.    && p > pmvi)  continue;
                if (p > 100e6 || (t > 1073.15 && p > 50e6))
                {
                        const auto
                    [ r, i ] = r6_inverse::density_pt (p, t, d, p * 1e-7, info::convergence);
                    if (i.converged)
                    {
                        d = r;
                        o   << t << " " << p << " " << d << "\n";
                        oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                        ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                        og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                        os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                        ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                        ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                        ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                        oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                        okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                        oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                        obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                        continue;
                    }
                    e << t << " " << p << " " << 1 << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                    const auto
                [ r, i ] = r6_inverse::density_pt (p, t, r7::density_pt (p, t), p * 1e-7, info::convergence);
                if (i.converged)
                {
                    d = r;
                    o   << t << " " << p << " " << d << "\n";
                    oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                    ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                    og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                    os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                    ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                    ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                    ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                    oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                    okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                    oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                    obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                    continue;
                }
                e << t << " " << p << " " << 1 << "\n";
                c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                for (auto&& [ v, f ,d ]: i.convergence)
                {
                    c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                }
                continue;
            }
        }
    }{
            auto
        o   = std::ofstream { "density_r6_inv_ylog.dat" };
            auto
        e   = std::ofstream { "density_r6_inv_ylog_failed.dat" };
            auto
        c   = std::ofstream { "density_r6_inv_ylog_failed_convergence.txt" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r6_inv_ylog.dat" };
            auto
        og  = std::ofstream { "massic_gibbs_free_energy_r6_inv_ylog.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r6_inv_ylog.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r6_inv_ylog.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r6_inv_ylog.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r6_inv_ylog.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r6_inv_ylog.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r6_inv_ylog.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r6_inv_ylog.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r6_inv_ylog.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r6_inv_ylog.dat" };
            const auto
        t0 = 173.15;
            const auto
        t1 = 1273;
            const auto
        dt = (t1 - t0) / (n - 1);
            const auto
        p0 = 1;
            const auto
        p1 = 1000e6;
            const auto
        dp = exp (log (p1 / p0) / (n - 1));
        for (auto it = 0; it < n; ++it)
        {
                const auto
            t = t0 + it * dt;
                const auto
            pmih = r14::ih::melting_pressure_t (t);
                const auto
            pmiii = r14::iii::melting_pressure_t (t);
                const auto
            pmv = r14::v::melting_pressure_t (t);
                const auto
            pmvi = r14::vi::melting_pressure_t (t);
                [[ maybe_unused ]]const auto
            pmvii = r14::vii::melting_pressure_t (t);
                const auto
            ps = r14::sublimation_pressure_t (t);
                auto
            d = 1000.;
            for (auto ip = 0; ip < n; ++ip)
            {
                    const auto
                p = p0 * pow (dp, ip);
                if (p <  611.657   && t <= 273.16 && p > ps) continue;
                if (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)  continue;
                if (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p > pmiii) continue;
                if (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p > pmv)   continue;
                if (p >= 632.4e6   && p <  2216e6    && t < 355.    && p > pmvi)  continue;
                if (p > 100e6 || (t > 1073.15 && p > 50e6))
                {
                        const auto
                    [r, i] = r6_inverse::density_pt (p, t, d, p * 1e-7, info::convergence);
                    if (i.converged)
                    {
                        d = r;
                        o   << t << " " << p << " " << d << "\n";
                        oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                        ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                        og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                        os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                        ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                        ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                        ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                        oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                        okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                        oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                        obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                        continue;
                    }
                    e << t << " " << p << " " << 1 << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                    const auto
                [ r, i ] = r6_inverse::density_pt (p, t, r7::density_pt (p, t), p * 1e-7, info::convergence);
                if (i.converged)
                {
                    d = r;
                    o   << t << " " << p << " " << d << "\n";
                    oh  << t << " " << p << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                    ou  << t << " " << p << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                    og  << t << " " << p << " " << r6::massic_gibbs_free_energy_dt (d, t)             << "\n";
                    os  << t << " " << p << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                    ocp << t << " " << p << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                    ocv << t << " " << p << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                    ow  << t << " " << p << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                    oav << t << " " << p << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                    okt << t << " " << p << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                    oap << t << " " << p << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                    obp << t << " " << p << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
                    continue;
                }
                e << t << " " << p << " " << 1 << "\n";
                c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                for (auto&& [ v, f ,d ]: i.convergence)
                {
                    c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                }
                continue;
            }
        }
    }
}
    auto
diagram_r7 (int n)
{
    {
            auto
        od  = std::ofstream { "density_r7_ylin.dat" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r7_ylin.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r7_ylin.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r7_ylin.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r7_ylin.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r7_ylin.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r7_ylin.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r7_ylin.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r7_ylin.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r7_ylin.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r7_ylin.dat" };
            const auto
        t0 = 273.15;
            const auto
        t1 = 2273.15;
            const auto
        dt = (t1 - t0) / n;
            const auto
        p0 = 0;
            const auto
        p1 = 100e6;
            const auto
        dp = (p1 - p0) / n;
        for (auto it = 0; it < n + 1; ++it)
        {
                const auto
            t = t0 + it * dt;
            for (auto ip = 1; ip < n + 1; ++ip)
            {
                    const auto
                p = p0 + ip * dp;
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                    od  << t << " " << p << " " << r7::density_pt (p, t)                              << "\n";
                    oh  << t << " " << p << " " << r7::massic_enthalpy_pt (p, t)                      << "\n";
                    ou  << t << " " << p << " " << r7::massic_internal_energy_pt (p, t)               << "\n";
                    os  << t << " " << p << " " << r7::massic_entropy_pt (p, t)                       << "\n";
                    ocp << t << " " << p << " " << r7::massic_isobaric_heat_capacity_pt (p, t)        << "\n";
                    ocv << t << " " << p << " " << r7::massic_isochoric_heat_capacity_pt (p, t)       << "\n";
                    ow  << t << " " << p << " " << r7::speed_of_sound_pt (p, t)                       << "\n";
                    oav << t << " " << p << " " << r7::isobaric_cubic_expansion_coefficient_pt (p, t) << "\n";
                    okt << t << " " << p << " " << r7::isothermal_compressibility_pt (p, t)           << "\n";
                    oap << t << " " << p << " " << r7::relative_pressure_coefficient_pt (p, t)        << "\n";
                    obp << t << " " << p << " " << r7::isothermal_stress_coefficient_pt (p, t)        << "\n";
                }
                catch (...) {}
            }
        }
    }
    {
            auto
        od  = std::ofstream { "density_r7_ylog.dat" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r7_ylog.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r7_ylog.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r7_ylog.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r7_ylog.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r7_ylog.dat" };
            auto
        ow = std::ofstream { "speed_of_sound_r7_ylog.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r7_ylog.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r7_ylog.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r7_ylog.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r7_ylog.dat" };
            const auto
        t0 = 273.15;
            const auto
        t1 = 2273.15;
            const auto
        dt = (t1 - t0) / n;
            const auto
        p0 = 1;
            const auto
        p1 = 100e6;
            const auto
        dp = exp (log (p1 / p0) / n);
        for (auto it = 0; it < n + 1; ++it)
        {
                const auto
            t = t0 + it * dt;
            for (auto ip = 1; ip < n + 1; ++ip)
            {
                    const auto
                p = p0 * pow (dp, ip);
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                    od  << t << " " << p << " " << r7::density_pt (p, t)                              << "\n";
                    oh  << t << " " << p << " " << r7::massic_enthalpy_pt (p, t)                      << "\n";
                    ou  << t << " " << p << " " << r7::massic_internal_energy_pt (p, t)               << "\n";
                    os  << t << " " << p << " " << r7::massic_entropy_pt (p, t)                       << "\n";
                    ocp << t << " " << p << " " << r7::massic_isobaric_heat_capacity_pt (p, t)        << "\n";
                    ocv << t << " " << p << " " << r7::massic_isochoric_heat_capacity_pt (p, t)       << "\n";
                    ow  << t << " " << p << " " << r7::speed_of_sound_pt (p, t)                       << "\n";
                    oav << t << " " << p << " " << r7::isobaric_cubic_expansion_coefficient_pt (p, t) << "\n";
                    okt << t << " " << p << " " << r7::isothermal_compressibility_pt (p, t)           << "\n";
                    oap << t << " " << p << " " << r7::relative_pressure_coefficient_pt (p, t)        << "\n";
                    obp << t << " " << p << " " << r7::isothermal_stress_coefficient_pt (p, t)        << "\n";
                }
                catch (...) {}
            }
        }
    }
}

    auto
diagram_r7_vs_r6_inv (int n)
{
    {
            auto
        ode = std::ofstream { "density_r7_vs_r6_inv_error_ylin.dat" };
            auto
        odi = std::ofstream { "density_r7_vs_r6_inv_iter_ylin.dat" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r7_vs_r6_inv_error_ylin.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r7_vs_r6_inv_error_ylin.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r7_vs_r6_inv_error_ylin.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylin.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylin.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r7_vs_r6_inv_error_ylin.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylin.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r7_vs_r6_inv_error_ylin.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r7_vs_r6_inv_error_ylin.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r7_vs_r6_inv_error_ylin.dat" };
            const auto
        t0 = 273.15;
            const auto
        t1 = 2273.15;
            const auto
        dt = (t1 - t0) / n;
            const auto
        p0 = 0;
            const auto
        p1 = 100e6;
            const auto
        dp = (p1 - p0) / n;
        for (auto it = 0; it < n + 1; ++it)
        {
                const auto
            t = t0 + it * dt;
            for (auto ip = 1; ip < n + 1; ++ip)
            {
                    const auto
                p = p0 + ip * dp;
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                        const auto
                    d7 = r7::density_pt (p, t);
                        const auto
                    [ d6, i ] = r6_inverse::density_pt (p, t, d7, p * 1e-8, info::iterations);
                        const auto
                    dd = fabs (d7 - d6);
                    ode << t << " " << p << " " << dd / d6 << "\n";
                    odi << t << " " << p << " " << i.iteration_count << "\n";
                        double a, b;
                    a = r6::massic_enthalpy_dt (d6, t);
                    b = r7::massic_enthalpy_pt (p, t);
                    oh  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_internal_energy_dt (d6, t);
                    b = r7::massic_internal_energy_pt (p, t);
                    ou  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_entropy_dt (d6, t);
                    b = r7::massic_entropy_pt (p, t);
                    os  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_isobaric_heat_capacity_dt (d6, t);
                    b = r7::massic_isobaric_heat_capacity_pt (p, t);
                    ocp << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_isochoric_heat_capacity_dt (d6, t);
                    b = r7::massic_isochoric_heat_capacity_pt (p, t);
                    ocv << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::speed_of_sound_dt (d6, t);
                    b = r7::speed_of_sound_pt (p, t);
                    ow  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isobaric_cubic_expansion_coefficient_dt (d6, t);
                    b = r7::isobaric_cubic_expansion_coefficient_pt (p, t);
                    oav << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isothermal_compressibility_dt (d6, t);
                    b = r7::isothermal_compressibility_pt (p, t);
                    okt << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::relative_pressure_coefficient_dt (d6, t);
                    b = r7::relative_pressure_coefficient_pt (p, t);
                    oap << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isothermal_stress_coefficient_dt (d6, t);
                    b = r7::isothermal_stress_coefficient_pt (p, t);
                    obp << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                }
                catch (...) {}
            }
        }
    }
    {
            auto
        oe = std::ofstream { "density_r7_vs_r6_inv_error_ylog.dat" };
            auto
        oi = std::ofstream { "density_r7_vs_r6_inv_iter_ylog.dat" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r7_vs_r6_inv_error_ylog.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r7_vs_r6_inv_error_ylog.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r7_vs_r6_inv_error_ylog.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylog.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylog.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r7_vs_r6_inv_error_ylog.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylog.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r7_vs_r6_inv_error_ylog.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r7_vs_r6_inv_error_ylog.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r7_vs_r6_inv_error_ylog.dat" };
            const auto
        t0 = 273.15;
            const auto
        t1 = 2273.15;
            const auto
        dt = (t1 - t0) / n;
            const auto
        p0 = 1;
            const auto
        p1 = 100e6;
            const auto
        dp = exp (log (p1 / p0) / n);
        for (auto it = 0; it < n + 1; ++it)
        {
                const auto
            t = t0 + it * dt;
            for (auto ip = 0; ip < n + 1; ++ip)
            {
                    const auto
                p = p0 * pow (dp, ip);
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                        const auto
                    d7 = r7::density_pt (p, t);
                        const auto
                    [ d6, i ] = r6_inverse::density_pt (p, t, d7, p * 1e-8, info::iterations);
                        const auto
                    dd = fabs (d7 - d6);
                    oe << t << " " << p << " " << dd / d6 << "\n";
                    oi << t << " " << p << " " << i.iteration_count << "\n";
                        double a, b;
                    a = r6::massic_enthalpy_dt (d6, t);
                    b = r7::massic_enthalpy_pt (p, t);
                    oh  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_internal_energy_dt (d6, t);
                    b = r7::massic_internal_energy_pt (p, t);
                    ou  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_entropy_dt (d6, t);
                    b = r7::massic_entropy_pt (p, t);
                    os  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_isobaric_heat_capacity_dt (d6, t);
                    b = r7::massic_isobaric_heat_capacity_pt (p, t);
                    ocp << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::massic_isochoric_heat_capacity_dt (d6, t);
                    b = r7::massic_isochoric_heat_capacity_pt (p, t);
                    ocv << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::speed_of_sound_dt (d6, t);
                    b = r7::speed_of_sound_pt (p, t);
                    ow  << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isobaric_cubic_expansion_coefficient_dt (d6, t);
                    b = r7::isobaric_cubic_expansion_coefficient_pt (p, t);
                    oav << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isothermal_compressibility_dt (d6, t);
                    b = r7::isothermal_compressibility_pt (p, t);
                    okt << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::relative_pressure_coefficient_dt (d6, t);
                    b = r7::relative_pressure_coefficient_pt (p, t);
                    oap << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                    a = r6::isothermal_stress_coefficient_dt (d6, t);
                    b = r7::isothermal_stress_coefficient_pt (p, t);
                    obp << t << " " << p << " " << fabs ((a - b) / a) << "\n";
                }
                catch (...) {}
            }
        }
    }

}
    auto
d_sat_g_t (double t)
    -> std::optional <double>
{
    if (t < 273.16 || t > 647.096) return {};
        const auto
    p = r7::saturation_pressure_t (t);
        const auto
    d7 = r7::r2::density_pt (p, t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, d7, p * 1e-7, info::convergence);
    if (!info.converged) return {};
    return d;
}
    auto
d_sat_l_t (double t)
    -> std::optional <double>
{
    if (t < 273.16 || t > 647.096) return {};
        const auto
    p = r7::saturation_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    return d;
}
    auto
d_sub_t (double t)
    -> std::optional <double>
{
    if (t < 200. || t > 273.16) return {};
        const auto
    p = r14::sublimation_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    if (d > 1.) return {};
    return d;
}
    auto
d_mel_ih_t (double t)
    -> std::optional <double>
{
    if (t < 251.165 || t > 273.16) return {};
        const auto
    p = r14::ih::melting_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    if (d < 1e3) return {};
    return d;
}
    auto
d_mel_iii_t (double t)
    -> std::optional <double>
{
    if (t < 251.165 || t > 256.164) return {};
        const auto
    p = r14::iii::melting_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    if (d < 0.) return {};
    return d;
}
    auto
d_mel_v_t (double t)
    -> std::optional <double>
{
    if (t < 256.164 || t > 273.31) return {};
        const auto
    p = r14::v::melting_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    if (d < 0. || d > 2.e3) return {};
    return d;
}
    auto
d_mel_vi_t (double t)
    -> std::optional <double>
{
    if (t < 273.31 || t > 355.) return {};
        const auto
    p = r14::vi::melting_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    return d;
}
    auto
d_mel_vii_t (double t)
    -> std::optional <double>
{
    if (t < 355. || t > 715.) return {};
        const auto
    p = r14::vii::melting_pressure_t (t);
        const auto
    [ d, info ] = r6_inverse::density_pt (p, t, 1000., p * 1e-7, info::convergence);
    if (!info.converged) return {};
    return d;
}
    auto
diagram_r6 (int n)
{
    {
            auto
        op  = std::ofstream { "pressure_r6_ylin.dat" };
            auto
        oh  = std::ofstream { "massic_enthalpy_r6_ylin.dat" };
            auto
        ou  = std::ofstream { "massic_internal_energy_r6_ylin.dat" };
            auto
        os  = std::ofstream { "massic_entropy_r6_ylin.dat" };
            auto
        ocp = std::ofstream { "massic_isobaric_heat_capacity_r6_ylin.dat" };
            auto
        ocv = std::ofstream { "massic_isochoric_heat_capacity_r6_ylin.dat" };
            auto
        ow  = std::ofstream { "speed_of_sound_r6_ylin.dat" };
            auto
        oav = std::ofstream { "isobaric_cubic_expansion_coefficient_r6_ylin.dat" };
            auto
        okt = std::ofstream { "isothermal_compressibility_r6_ylin.dat" };
            auto
        oap = std::ofstream { "relative_pressure_coefficient_r6_ylin.dat" };
            auto
        obp = std::ofstream { "isothermal_stress_coefficient_r6_ylin.dat" };
            const auto
        t0 = 273.16;
            const auto
        t1 = 1273.16;
            const auto
        dt = (t1 - t0) / (n - 1);
            const auto
        d0 = 0.;
            const auto
        d1 = 1000.;
            const auto
        dd = (d1 - d0) / (n - 1);
        for (auto it = 0; it < n; ++it)
        {
                const auto
            t = t0 + it * dt;
                const auto
            d_sat_g = d_sat_g_t (t);
                const auto
            d_sat_l = d_sat_l_t (t);
                const auto
            d_sub = d_sub_t (t);
                const auto
            d_mel_ih = d_mel_ih_t (t);
                const auto
            d_mel_iii = d_mel_iii_t (t);
                const auto
            d_mel_v = d_mel_v_t (t);
                const auto
            d_mel_vi = d_mel_vi_t (t);
                const auto
            d_mel_vii = d_mel_vii_t (t);
            for (auto id = 1; id < n; ++id)
            {
                    const auto
                d = d0 + id * dd;
                    const auto
                p = r6::pressure_dt (d, t);
                if (d_mel_vii && d > *d_mel_vii) continue;
                if (d_mel_vi  && d > *d_mel_vi)  continue;
                op  << d << " " << t << " " << p                                                  << "\n";
                oh  << d << " " << t << " " << r6::massic_enthalpy_dt (d, t)                      << "\n";
                ou  << d << " " << t << " " << r6::massic_internal_energy_dt (d, t)               << "\n";
                os  << d << " " << t << " " << r6::massic_entropy_dt (d, t)                       << "\n";
                ocp << d << " " << t << " " << r6::massic_isobaric_heat_capacity_dt (d, t)        << "\n";
                ocv << d << " " << t << " " << r6::massic_isochoric_heat_capacity_dt (d, t)       << "\n";
                ow  << d << " " << t << " " << r6::speed_of_sound_dt (d, t)                       << "\n";
                oav << d << " " << t << " " << r6::isobaric_cubic_expansion_coefficient_dt (d, t) << "\n";
                okt << d << " " << t << " " << r6::isothermal_compressibility_dt (d, t)           << "\n";
                oap << d << " " << t << " " << r6::relative_pressure_coefficient_dt (d, t)        << "\n";
                obp << d << " " << t << " " << r6::isothermal_stress_coefficient_dt (d, t)        << "\n";
            }
        }
    }
}
    template <class F>
    auto
plot_dt_6 (
      std::string const& name
    , F&& f
    , double t0
    , double t1
    , int n
    , double init
    , std::optional <double> low = {}
    , std::optional <double> high = {}
){
    if (!low)  low  = std::numeric_limits <double>::lowest ();
    if (!high) high = std::numeric_limits <double>::max ();
        const auto
    dt = (t1 - t0) / (n - 1);
        auto
    map = std::map <double, std::pair <double, double>> {};
    for (auto i = 0; i != n; ++i)
    {
            const auto
        t = t0 + dt * i;
            const auto
        p = std::forward <F> (f) (t);
            const auto
        [ d, info ] = r6_inverse::density_pt (p, t, /*1000.*/init, p * 1e-7, info::convergence);
        if (!info.converged) continue;
        if (d < low || d > high) continue;
        map.insert ({ d, { t, p } });
    }
        auto
    o = std::ofstream { name };
    for (auto&& [d, v]: map)
    {
            auto&&
        [ t, p ] = v;
        o << d << " " << t << " " << p << "\n";
    }
}
    template <class F>
    auto
plot_dt_6_7 (
      std::string const& name
    , F&& f
    , double t0
    , double t1
    , int n
    , std::optional <double> low = {}
    , std::optional <double> high = {}
){
    if (!low)  low  = std::numeric_limits <double>::lowest ();
    if (!high) high = std::numeric_limits <double>::max ();
        const auto
    dt = (t1 - t0) / (n - 1);
        auto
    map = std::map <double, std::pair <double, double>> {};
    for (auto i = 0; i != n; ++i)
    {
            const auto
        t = t0 + dt * i;
            const auto
        p = std::forward <F> (f) (t);
            const auto
        d7 = r7::r2::density_pt (p, t);
            const auto
        [ d, info ] = r6_inverse::density_pt (p, t, d7, p * 1e-7, info::convergence);
        if (!info.converged) continue;
        if (d < low || d > high) continue;
        map.insert ({ d, { t, p } });
    }
        auto
    o = std::ofstream { name };
    for (auto&& [d, v]: map)
    {
            auto&&
        [ t, p ] = v;
        o << d << " " << t << " " << p << "\n";
    }
}
    auto
lines_dt (int n)
{
    plot_dt_6_7 ("saturation_gaz_dt.dat", r7::saturation_pressure_t    <double>, 273.16, 647.096, n);
    plot_dt_6   ("saturation_liq_dt.dat", r7::saturation_pressure_t    <double>, 273.16, 647.096, n, 1000.);
    plot_dt_6   ("sublimation_dt.dat",    r14::sublimation_pressure_t  <double>, 200., 273.16, n, 1., {}, 1.);
    plot_dt_6   ("melting_ih_dt.dat",     r14::ih::melting_pressure_t  <double>, 251.165, 273.16,  n, 1000., 1e3);
    plot_dt_6   ("melting_iii_dt.dat",    r14::iii::melting_pressure_t <double>, 251.165, 256.164, n, 1000., 0.);
    plot_dt_6   ("melting_v_dt.dat",      r14::v::melting_pressure_t   <double>, 256.164, 273.31,  n, 1000., 0., 2e3);
    plot_dt_6   ("melting_vi_dt.dat",     r14::vi::melting_pressure_t  <double>, 273.31, 355., n, 1000.);
    plot_dt_6   ("melting_vii_dt.dat",    r14::vii::melting_pressure_t <double>, 355.,   715., n, 1000.);
}
    template <class F>
    auto
plot_pt (std::string const& name, F&& f, double t0, double t1, int n = 100)
{
        auto
    o = std::ofstream { name };
        const auto
    dt = (t1 - t0) / (n);
    for (auto i = 0; i != n + 1; ++i)
    {
            const auto
        t = t0 + i * dt;
        o << t << " " << std::forward <F> (f) (t) << "\n";
    }
}
    auto
lines_pt (int n)
{
    plot_pt ("saturation.dat",  r7::saturation_pressure_t    <double>, 273.16,  647.096, n);
    plot_pt ("sublimation.dat", r14::sublimation_pressure_t  <double>, 200.,    273.16,  n);
    plot_pt ("melting_ih.dat",  r14::ih::melting_pressure_t  <double>, 251.165, 273.16,  n);
    plot_pt ("melting_iii.dat", r14::iii::melting_pressure_t <double>, 251.165, 256.164, n);
    plot_pt ("melting_v.dat",   r14::v::melting_pressure_t   <double>, 256.164, 273.31,  n);
    plot_pt ("melting_vi.dat",  r14::vi::melting_pressure_t  <double>, 273.31,  355.,    n);
    plot_pt ("melting_vii.dat", r14::vii::melting_pressure_t <double>, 355.,    715.,    n);
}
    int
main ()
{
        constexpr auto
    n = 501;
        auto
    work = std::vector <std::thread> {};
    work.emplace_back (lines_pt, n);
    work.emplace_back (lines_dt, n);
    work.emplace_back (diagram_r7, n);
    work.emplace_back (diagram_r6_inv, n);
    work.emplace_back (diagram_r6_inv_ext, n);
    work.emplace_back (diagram_r7_vs_r6_inv, n);
    work.emplace_back (diagram_r6, n);
    for (auto&& t: work)
    {
        t.join ();
    }
    return 0;
}
