#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r6_inverse.hpp"
#include "../include/isto/iapws/r7.hpp"
#include "../include/isto/iapws/r10.hpp"
#include "../include/isto/iapws/r14.hpp"
#include "../include/isto/iapws/g12.hpp"
    using namespace isto::iapws;
#include <fstream>
#include <fmt/format.h>
    using fmt::print, fmt::format;
#include <optional>
#include <map>

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
        [ d, info ] = r6_inverse::density_pt (p, t, /*1000.*/init, info::convergence);
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
){
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
        [ d, info ] = r6_inverse::density_pt (p, t, d7, info::convergence);
        if (!info.converged) continue;
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
plot_dt_6_p (
      std::string const& name
    , F&& f
    , double p0
    , double p1
    , int n
    , double init
    , std::optional <double> low = {}
    , std::optional <double> high = {}
){
    if (!low)  low  = std::numeric_limits <double>::lowest ();
    if (!high) high = std::numeric_limits <double>::max ();
        const auto
    dp = (p1 - p0) / (n - 1);
        auto
    map = std::map <double, std::pair <double, double>> {};
    for (auto i = 0; i != n; ++i)
    {
            const auto
        p = p0 + dp * i;
            const auto
        t = std::forward <F> (f) (p);
            const auto
        d = r10::density_pt (p, t);
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
    plot_dt_6_7 ("saturation_gaz_dt.dat", r7::saturation_pressure_t    <double>, 273.15,  647.096, n);
    plot_dt_6   ("saturation_liq_dt.dat", r7::saturation_pressure_t    <double>, 273.15,  647.096, n, 1000.);
    plot_dt_6   ("sublimation_dt.dat",    r14::sublimation_pressure_t  <double>, 50.,     273.16,  n, 1., {}, 1.);
    plot_dt_6   ("melting_ih_dt.dat",     r14::ih::melting_pressure_t  <double>, 251.165, 273.16,  n, 1000., 1e3);
    plot_dt_6   ("melting_iii_dt.dat",    r14::iii::melting_pressure_t <double>, 251.165, 256.164, n, 1000., 0.);
    plot_dt_6   ("melting_v_dt.dat",      r14::v::melting_pressure_t   <double>, 256.164, 273.31,  n, 1000., 0., 2e3);
    plot_dt_6   ("melting_vi_dt.dat",     r14::vi::melting_pressure_t  <double>, 273.31,  355.,    n, 1000.);
    plot_dt_6   ("melting_vii_dt.dat",    r14::vii::melting_pressure_t <double>, 355.,    715.,    n, 1000.);
    plot_dt_6_p ("ice_nucleation_dt.dat", g12::homogeneous_ice_nucleation_limit_temperature_p <double>, 1., 1500e6, n, 1000.);
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
    template <class F>
    auto
plot_pt_p (std::string const& name, F&& f, double p0, double p1, int n = 100)
{
        auto
    o = std::ofstream { name };
        const auto
    dp = exp (log (p1 / p0) / n);
    for (auto i = 0; i != n + 1; ++i)
    {
            const auto
        p = p0 + pow (dp, i);
        o << std::forward <F> (f) (p) << " " << p << "\n";
    }
}
    auto
lines_tp (int n)
{
    plot_pt   ("saturation_tp.dat",     r7::saturation_pressure_t    <double>, 273.15,  647.096, n);
    plot_pt   ("sublimation_tp.dat",    r14::sublimation_pressure_t  <double>, 50.,     273.16,  n);
    plot_pt   ("melting_ih_tp.dat",     r14::ih::melting_pressure_t  <double>, 251.165, 273.16,  n);
    plot_pt   ("melting_iii_tp.dat",    r14::iii::melting_pressure_t <double>, 251.165, 256.164, n);
    plot_pt   ("melting_v_tp.dat",      r14::v::melting_pressure_t   <double>, 256.164, 273.31,  n);
    plot_pt   ("melting_vi_tp.dat",     r14::vi::melting_pressure_t  <double>, 273.31,  355.,    n);
    plot_pt   ("melting_vii_tp.dat",    r14::vii::melting_pressure_t <double>, 355.,    715.,    n);
    plot_pt_p ("ice_nucleation_tp.dat", g12::homogeneous_ice_nucleation_limit_temperature_p <double>, 1., 1500e6, n);
}
    auto
make_gpl ()
{
    {
            auto
        f = std::ofstream { "lines_tp.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960 enhanced
set xlabel "temperature [K]"
set ylabel "pressure [MPa]"
set cblabel "density [kg/m3]"
set ytics format "%.0e"
set cbtics format "%.0e"
set logscale y
set output "lines_tp.png"
set key bottom
set title "Phases boundaries in the temperature-pressure plane"
set xrange [173:773]
set yrange [1e-6:1e6]

set object circle center 251.165,208.556 size 3
set label "Triple point L-Ih-III: 251.165 K, 208.566 MPa" at 280,150 left
set arrow from 280,150 to 255,208

set object circle center 256.164,350.1 size 3
set label "Triple point L-III-V: 256.164 K, 350.1 MPa" at 300,400 left
set arrow from 300,400 to 260,350.1

set object circle center 273.31,632.4 size 3
set label "Triple point L-V-VI: 273.31 K, 632.4 MPa" at 270,2e4 center
set arrow from 270,1.5e4 to 273,800

set object circle center 355,2216 size 3
set label "Triple point L-VI-VII: 355 K, 2216 MPa" at 355,1300 left

set object circle center 647.096,22.064 size 3
set label "Critical point: 647.096 K, 22.064 MPa" at 647.096,40 center

set object circle center 273.16,611.657e-6 size 3
set label "Triple point S-L-G: 273.16 K, 611.657 Pa" at 280,500e-6 left

plot \
      "sublimation_tp.dat"    using 1:($2*1e-6) w l t "sublimation line" \
    , "saturation_tp.dat"     using 1:($2*1e-6) w l t "saturation line"  \
    , "melting_ih_tp.dat"     using 1:($2*1e-6) w l t "melting line Ih"  \
    , "melting_iii_tp.dat"    using 1:($2*1e-6) w l t "melting line III" \
    , "melting_v_tp.dat"      using 1:($2*1e-6) w l t "melting line V"   \
    , "melting_vi_tp.dat"     using 1:($2*1e-6) w l t "melting line VI"  \
    , "melting_vii_tp.dat"    using 1:($2*1e-6) w l t "melting line VII" \
    , "ice_nucleation_tp.dat" using 1:($2*1e-6) w l t "homogeneous ice nucleation line"


# vim:ft=gnuplot
    )";
    }
    {
            auto
        f = std::ofstream { "lines_dt.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960 enhanced
set ylabel "temperature [K]"
set xlabel "density [kg/m3]"
set output "lines_dt.png"
set title "Phases boundaries in the density-temperature plane"
set xrange [-100:1900]
set yrange [150:700]
set key top center

set object circle center 322,647.096 size 3
set label "Critical point: 322 kg/m^3, 647.096 K" at 322,670 center

set arrow from -100,273.16 to 1900,273.16 nohead dt 2
set style textbox opaque fillcolor "white" noborder
set label "Triple point temperature 273.16 K" at 500,273.16 center front boxed

plot                                                        \
      "sublimation_dt.dat"    w l t "Sublimation line"      \
    , "saturation_gaz_dt.dat" w l t 'Saturation line (gaz)' \
    , "saturation_liq_dt.dat" w l t 'Saturation line (liq)' \
    , "melting_ih_dt.dat"     w l t "Melting line Ih"       \
    , "melting_iii_dt.dat"    w l t "Melting line III"      \
    , "melting_v_dt.dat"      w l t "Melting line V"        \
    , "melting_vi_dt.dat"     w l t "Melting line VI"       \
    , "melting_vii_dt.dat"    w l t "Melting line VII"      \
    , "ice_nucleation_dt.dat" w l t "homogeneous ice nucleation line"

# vim:ft=gnuplot
    )";
    }
}

    int
main ()
{
        constexpr auto
    n = 501;
    lines_tp (n);
    lines_dt (n);
    make_gpl ();
    return 0;
}
