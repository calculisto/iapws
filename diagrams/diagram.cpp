#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r7.hpp"
#include "../include/isto/iapws/r10.hpp"
#include "../include/isto/iapws/r14.hpp"
    using namespace isto::iapws;
#include <fstream>
#include <fmt/format.h>
    using fmt::print;

    template <class F>
    auto
plot (std::string const& name, F&& f, double t0, double t1, int n = 100)
{
        auto
    o = std::ofstream { name };
        auto const
    dt = (t1 - t0) / (n);
    for (auto i = 0; i != n + 1; ++i)
    {
            auto const
        t = t0 + i * dt;
        o << t << " " << std::forward <F> (f) (t) << "\n";
    }
}

    auto
diagram_r6 (int n = 500)
{
    {
            auto
        o = std::ofstream { "density_r6.dat" };
            auto const
        t0 = 273.15;
            auto const
        t1 = 1273;
            auto const
        dt = (t1 - t0) / n;
            auto const
        p0 = 0;
            auto const
        p1 = 1000e6;
            auto const
        dp = (p1 - p0) / n;
        for (auto it = 0; it < n + 1; ++it)
        {
                auto const
            t = t0 + it * dt;
            for (auto ip = 0; ip < n + 1; ++ip)
            {
                    auto const
                p = p0 + ip * dp;
                if (t > 1073.15 && p > 50e6) continue;
                    auto const
                d = r6_inverse::density_pt (p, t);
                o << t << " " << p << " " << d << "\n";
                catch (...) {}
            }
        }
    }
}
    auto
diagram_r7 (int n = 500)
{
    {
            auto
        o = std::ofstream { "density_r7.dat" };
            auto const
        t0 = 273.15;
            auto const
        t1 = 2273.15;
            auto const
        dt = (t1 - t0) / n;
            auto const
        p0 = 0;
            auto const
        p1 = 100e6;
            auto const
        dp = (p1 - p0) / n;
        for (auto it = 0; it < n + 1; ++it)
        {
                auto const
            t = t0 + it * dt;
            for (auto ip = 0; ip < n + 1; ++ip)
            {
                    auto const
                p = p0 + ip * dp;
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                        auto const
                    d = r7::density_pt (p, t);
                    o << t << " " << p << " " << d << "\n";
                }
                catch (...) {}
            }
        }
    }
    {
            auto
        o = std::ofstream { "density_r7_plog.dat" };
            auto const
        t0 = 273.15;
            auto const
        t1 = 2273.15;
            auto const
        dt = (t1 - t0) / n;
            auto const
        p0 = 1;
            auto const
        p1 = 100e6;
            auto const
        dp = exp (log (p1 / p0) / n);
        for (auto it = 0; it < n + 1; ++it)
        {
                auto const
            t = t0 + it * dt;
            for (auto ip = 0; ip < n + 1; ++ip)
            {
                    auto const
                p = p0 * pow (dp, ip);
                if (t > 1073.15 && p > 50e6) continue;
                try
                {
                        auto const
                    d = r7::density_pt (p, t);
                    o << t << " " << p << " " << d << "\n";
                }
                catch (...) {}
            }
        }
    }
}

    auto
diagram (int n = 500)
{
    {
            auto const
        t0 = 173.15;
            auto const
        t1 = 1073.15;
            auto const
        dt = (t1 - t0) / (n);
            auto const
        p0 = 1;
            auto const
        p1 = 1e8;
            auto const
        dp = exp (log (p1 / p0) / n);
            auto
        o = std::ofstream { "density.dat" };
        for (auto it = 0; it < n + 1; ++it)
        {
                auto const
            t = t0 + it * dt;
                auto const
            pm = r14::ih::melting_pressure_t (t);
                auto const
            ps = r14::sublimation_pressure_t (t);
            for (auto ip = 0; ip < n + 1; ++ip)
            {
                    auto const
                p = p0 * pow (dp, ip);
                if (p < pm && p > ps)
                {
                        auto const
                    d = r10::density_pt (p, t);
                    o << t << " " << p << " " << d << "\n";
                    continue;
                }
                    auto const
                d = r7::density_pt (p, t);
                o << t << " " << p << " " << d << "\n";
            }
        }
    }
}
    int
main ()
{
    diagram_r7 ();
    plot ("saturation.dat",  r7::saturation_pressure_t <double>, 273.16, 647.096);
    plot ("sublimation.dat", r14::sublimation_pressure_t <double>, 200., 273.16);
    plot ("melting_ih.dat",  r14::ih::melting_pressure_t <double>, 251.165, 273.16);
    plot ("melting_iii.dat", r14::iii::melting_pressure_t <double>, 251.165, 256.164);
    plot ("melting_v.dat",   r14::v::melting_pressure_t <double>, 256.164, 273.31);
    plot ("melting_vi.dat",  r14::vi::melting_pressure_t <double>, 273.31, 355.);
    plot ("melting_vii.dat", r14::vii::melting_pressure_t <double>, 355., 715.);
    //diagram ();
    return 0;
}
