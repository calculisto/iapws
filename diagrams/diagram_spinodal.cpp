#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6.hpp"
    using namespace isto::iapws;
#include <isto/root_finding/root_finding.hpp>
    using namespace isto::root_finding;
#include <fmt/format.h>
    using fmt::print, fmt::format;
#include <fstream>
#include "../include/isto/iapws/r7.hpp" 
#include "../include/isto/iapws/r6_inverse.hpp" 

    auto
saturation_density_gaz_t (double temperature)
{
        const auto
    pressure = r7::saturation_pressure_t (temperature);
        const auto
    density_r7 = r7::r2::density_pt (pressure, temperature);
        const auto
    density = r6_inverse::density_pt (
          pressure
        , temperature
        , density_r7
        , pressure * 1e-7
    );
    return density;
}

    auto
f (double delta, double tau)
{
        using namespace r6::r6_95_2016::detail;
    return 1. + 2. * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau);
}
    auto
df (double delta, double tau)
{
        using namespace r6::r6_95_2016::detail;
    return (delta < 1. ? std::nan ("") : 2. * phi_r_d (delta, tau) + 4. * delta * phi_r_dd (delta, tau) + delta * delta * phi_r_ddd (delta, tau));
}
    auto
g (double delta, double tau)
{
        using namespace r6_gas::r6_95_2016::detail;
    return 1. + 2. * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau);
}
    auto
dg (double delta, double tau)
{
        using namespace r6_gas::r6_95_2016::detail;
    return 2. * phi_r_d (delta, tau) + 4. * delta * phi_r_dd (delta, tau) + delta * delta * phi_r_ddd (delta, tau);
}

    auto
plot (double tau, int index)
{
        auto
    o = std::ofstream { format ("spi_{}.dat", index) };
        const auto
    d0 = 0.;
        const auto
    d1 = 2;
        const auto
    n = 1001;
        const auto
    dd = (d1 - d0) / (n - 1);
    for (auto i = 0; i < n; ++i)
    {
            const auto
        delta = d0 + i * dd;
        o 
            << delta 
            << " " << f  (delta, tau) 
            << " " << df (delta, tau)
            << " " << g  (delta, tau)
            << " " << dg (delta, tau)
        << "\n";
    }
}

    auto
make_gpl ()
{
    {
            auto
        f = std::ofstream { "spinodal_tp.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960 enhanced
set xlabel "temperature [K]"
set ylabel "pressure [MPa]"
set cblabel "density [kg/m3]"
set ytics format "%.0e"
set cbtics format "%.0e"
set output "spinodal_tp.png"
set key bottom
set title "Phases boundaries in the temperature-pressure plane"
set xrange [173:873]
set yrange [-200:40]

set object circle center 647.096,22.064 size 3
set label "Critical point: 647.096 K, 22.064 MPa" at 647.096,30 center

set object circle center 273.16,611.657e-6 size 3
set label "Triple point S-L-G: 273.16 K, 611.657 Pa" at 273,-8 left

plot \
      "spinodal_liq.dat"   using 2:($3*1e-6) w l t 'Spinodal line (liq)'  \
    , "spinodal_gaz_candidates.dat" u 1:($3*1e-6) w l t"Spinodal line (gas), candidate 1" \
    , "spinodal_gaz_candidates.dat" u 1:($5*1e-6) w l t"Spinodal line (gas), candidate 2" \
    , "spinodal_gaz.dat"      u 2:($3*1e-6) w l t "Spinodal line (gaz), R6 gas" \
    , "sublimation_tp.dat" using 1:($2*1e-6) w l t "sublimation line" dt 2 \
    , "saturation_tp.dat"  using 1:($2*1e-6) w l t "saturation line"  dt 2 \
    , "melting_ih_tp.dat"  using 1:($2*1e-6) w l t "melting line Ih"  dt 2 \

# vim:ft=gnuplot
    )";
    }
    {
            auto
        f = std::ofstream { "spinodal_dt.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960 enhanced
set ylabel "temperature [K]"
set xlabel "density [kg/m3]"
set output "spinodal_dt.png"
set key center center
set title "Phases boundaries in the density-temperature plane"
set xrange [-100:1100]
set yrange [200:700]

set object circle center 322,647.096 size 3
set label "Critical point: 322 kg/m^3, 647.096 K" at 322,670 center

set arrow from -100,273.16 to 1900,273.16 nohead dt 2
set style textbox opaque fillcolor "white" noborder
set label "Triple point temperature 273.16 K" at 500,273.16 center front boxed

plot                                                                              \
      "spinodal_liq.dat"      w l t 'Spinodal line (liq)'                         \
    , "spinodal_gaz_candidates.dat" u 2:1 w l t"Spinodal line (gaz), candidate 1" \
    , "spinodal_gaz_candidates.dat" u 4:1 w l t"Spinodal line (gaz), candidate 2" \
    , "spinodal_gaz.dat"      w l t "Spinodal line (gas), R6 gas"                 \
    , "sublimation_dt.dat"    w l t "Sublimation line"      dt 2                  \
    , "saturation_gaz_dt.dat" w l t 'Saturation line (gas)' dt 2                  \
    , "saturation_liq_dt.dat" w l t 'Saturation line (liq)' dt 2                  \

# vim:ft=gnuplot
    )";
    }
}

    int
main ()
{
    {
            auto
        o = std::ofstream { "spinodal_liq.dat" };
            const auto
        T0 = 647.096;
            const auto
        T1 = 173.16;
            const auto
        n = 1001;
            const auto
        dT = (T1 - T0) / (n - 1);
            auto
        d_init = 1.1;
        for (auto it = 0; it != n; ++it)
        {
                const auto
            T = T0 + it * dT;
                const auto
            tau = 647.096 / T;
                const auto
            [ delta, info ] = newton (
                  [=](auto delta){ return f (delta, tau); }
                , [=](auto delta){ return df (delta, tau); }
                , d_init
                , [=](auto x)
                  { 
                        using std::abs;
                    return abs (x) < 1e-10; 
                  }
                , {} // options
                , info::convergence
            );
            if (!info.converged)
            {
                print ("Failed at tau= {}\n", tau);
                plot (tau, it);
                for (auto [v, f, df]: info.convergence)
                {
                    print ("{}, {}, {}\n", v, f, df);
                }
                continue;
            }
                const auto
            D = delta * 322.;
                const auto
            P = r6::pressure_dt (D, T);
            o << D << " " << T << " " << P << "\n";
            d_init = delta;
        }
    }
    {
            auto
        o = std::ofstream { "spinodal_gaz.dat" };
            auto
        p = std::ofstream { "spinodal_gaz_candidates.dat" };
            const auto
        T0 = 647.096;
            const auto
        T1 = 173.16;
            const auto
        n = 101;
            const auto
        dT = (T1 - T0) / (n - 1);
        p << 647.096 << " " << 322. << " " << 22.064e6 << " " << 322. << " " << 22.064e6 << "\n";
        for (auto it = 1; it != n; ++it)
        {
                const auto
            T = T0 + it * dT;
                const auto
            tau = 647.096 / T;
                const auto
            saturation_density = saturation_density_gaz_t (T);
                const auto
            sd = saturation_density / 322.;
            // ===
            {
                    const auto
                [ delta, info ] = newton (
                      [=](auto delta){ return g (delta, tau); }
                    , [=](auto delta){ return dg (delta, tau); }
                    , sd
                    , [=](auto x)
                      { 
                            using std::abs;
                        return abs (x) < 1e-10; 
                      }
                    , {} // options
                    , info::convergence
                );
                if (!info.converged)
                {
                    /*
                    print ("Gas failed at tau= {}\n", tau);
                    plot (tau, it);
                    for (auto [v, f, df]: info.convergence)
                    {
                        print ("{}, {}, {}\n", v, f, df);
                    }
                    continue;
                    */
                }
                if (info.converged)
                {
                        const auto
                    D = delta * 322.;
                        const auto
                    P = r6_gas::pressure_dt (D, T);
                    o << D << " " << T << " " << P << "\n";
                }
                //d_init = delta;
            }
            // ===
            // ---
                const auto
            dd = (1. - sd) / 100;
                auto
            d0 = sd;
                auto
            f0 = f (d0, tau);
            p << 647.096 / tau;
            if (f0 == 0.) p << " " << d0 * 322.;
            for (auto j = 0; j < 100; ++j)
            {
                    const auto
                d1 = d0 + dd;
                    const auto
                f1 = f(d1, tau);
                if (f1 == 0.) p << " " << d1 * 322.;
                if (f0 * f1 < 0.)
                {
                        const auto
                    [ delta, info ] = zhang (
                          [=](auto delta){ return f (delta, tau); }
                        , d0
                        , d1
                        , [=](auto a, auto b, auto fa, auto fb)
                          { 
                                    constexpr double
                                tol = 1e-10;
                                return 
                                       fa == 0.0
                                    || fb == 0.0
                                    || fabs (b - a) < tol
                                ;
                          }
                        , {} // options
                        , info::convergence
                    );
                    if (!info.converged)
                    {
                        print ("Failed at tau= {} ({}) brackets are  [{}, {}]\n", tau, it, d0, d1);
                        plot (tau, it);
                        for (auto [a, b, fa, fb]: info.convergence)
                        {
                            print ("{}, {}, {}, {}\n", a, b, fa, fb);
                        }
                        continue;
                    }
                        const auto
                    D = delta * 322.;
                        const auto
                    P = r6::pressure_dt (D, T);
                    p << " " << D << " " << P;
                }
                d0 = d1;
                f0 = f1;
            }
            p << "\n";
            // ---

            /*
                const auto
            w = 1. - sd;
                const auto
            //da = sd + 0.0 * w;
            da = sd + 0.5 * w;
                const auto
            //db = da + 0.5 * w;
            db = sd + 0.9 * w;
                const auto
            [ delta, info ] = zhang (
                  [=](auto delta){ return 1. + 2. * delta * phi_r_d (delta, tau) + delta * delta * phi_r_dd (delta, tau); }
                , da
                , db
                , [=](auto a, auto b, auto fa, auto fb)
                  { 
                            constexpr double
                        tol = 1e-10;
                        return 
                               fa == 0.0
                            || fb == 0.0
                            || fabs (b - a) < tol
                        ;
                  }
                , {} // options
                , info::convergence
            );
            if (!info.converged)
            {
                print ("Failed at tau= {} ({}) brackets are  [{}, {}]\n", tau, it, da, db);
                plot (tau, it);
                for (auto [a, b, fa, fb]: info.convergence)
                {
                    print ("{}, {}, {}, {}\n", a, b, fa, fb);
                }
                continue;
            }
                const auto
            D = delta * 322.;
            o << D << " " << T << "\n";
            */
        }
    }
    make_gpl ();
}
