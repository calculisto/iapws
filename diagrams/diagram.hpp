#pragma once
#define ISTO_IAPWS_FORCE_RELAXED 1
#include <functional>
#include <string>
#include <unordered_map>
#include <fstream>
#include <fmt/format.h>
    using fmt::print, fmt::format;
    using namespace std::literals;
#include "../include/isto/iapws/r7.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
#include "../include/isto/iapws/r14.hpp"
#include "../include/isto/iapws/g12.hpp"
    using namespace isto::iapws;

    struct
quantity_t
{
        std::string
    label;
        std::string
    unit;
        double
    scale;
};
    auto const
quantities = std::unordered_map 
{
       std::pair
      { "density"s,                        quantity_t { "Density"s,                                                              "kg/m^3"s,   1.,  }}
    , { "pressure",                                   { "Pressure",                                                              "MPa",       1e-6,}}
    , { "temperature",                                { "Temperature",                                                           "K",         1.,  }}
    , { "massic_enthalpy",                            { "Massic enthalpy",                                                       "kJ/kg",     1e3, }}
    , { "massic_entropy",                             { "Massic entropy",                                                        "kJ/kg/K",   1e3, }}
    , { "massic_internal_energy",                     { "Massic internal energy",                                                "kJ/kg",     1e3, }}
    , { "massic_gibbs_free_energy",                   { "Massic Gibbs free energy",                                              "kJ/kg",     1e3, }}
    , { "massic_isobaric_heat_capacity",              { "Massic isobaric heat capacity",                                         "kJ/kg/K",   1e3, }}
    , { "massic_isochoric_heat_capacity",             { "Massic isochoric heat capacity",                                        "kJ/kg/K",   1e3, }}
    , { "speed_of_sound",                             { "Speed of sound",                                                        "m/s",       1.,  }}
    , { "isobaric_cubic_expansion_coefficient",       { "Isobaric cubic expansion coefficient",                                  "K^{ -1}",   1.,  }}
    , { "thermal_expansion_coefficient",              { "Thermal expansion coefficient",                                         "K^{ -1}",   1.,  }}
    , { "isothermal_compressibility",                 { "Isothermal compressibility",                                            "MPa^{ -1}", 1e6, }}
    , { "isothermal_stress_coefficient",              { "Isothermal stress coefficient",                                         "kg/m^3",    1.,  }}
    , { "relative_pressure_coefficient",              { "Relative pressure coefficient",                                         "K^{ -1}",   1.,  }}
    , { "relative_pressure_coefficient",              { "Relative pressure coefficient",                                         "K^{ -1}",   1.,  }}
    , { "iteration_count",                            { "Iteration count",                                                       "",          1.,  }}
    , { "delta_density",                              { "Density difference, relative to R6 value",                              "",          1.,  }}
    , { "delta_massic_enthalpy",                      { "Massic enthalpy difference, relative to R6 value",                      "",          1.,  }}
    , { "delta_massic_entropy",                       { "Massic entropy difference, relative to R6 value",                       "",          1.,  }}
    , { "delta_massic_internal_energy",               { "Massic internal energy difference, relative to R6 value",               "",          1.,  }}
    , { "delta_massic_isobaric_heat_capacity",        { "Massic isobaric heat capacity difference, relative to R6 value",        "",          1.,  }}
    , { "delta_massic_isochoric_heat_capacity",       { "Massic isochoric heat capacity difference, relative to R6 value",       "",          1.,  }}
    , { "delta_speed_of_sound",                       { "Speed of sound difference, relative to R6 value",                       "",          1.,  }}
    , { "delta_isobaric_cubic_expansion_coefficient", { "Isobaric cubic expansion coefficient difference, relative to R6 value", "",          1.,  }}
    , { "delta_isothermal_compressibility",           { "Isothermal compressibility difference, relative to R6 value",           "",          1.,  }}
    , { "delta_isothermal_stress_coefficient",        { "Isothermal stress coefficient difference, relative to R6 value",        "",          1.,  }}
    , { "delta_relative_pressure_coefficient",        { "Relative pressure coefficient difference, relative to R6 value",        "",          1.,  }}
    , { "viscosity",                                  { "Viscosity",                                                             "Pa.s",      1.,  }}
};

    class
exclusion_base_t
{
public:
        virtual void
    init (double)
    {}
        virtual bool
    operator () (double, double) const
    {
        return false;
    }
        virtual
    ~exclusion_base_t () = default;
};
    class
exclusion_dt_t
    : public exclusion_base_t
{
        auto
    d_sat_g_t (double t)
    {
        if (t < 273.16 || t > 647.096) return 0.;
            const auto
        p = r7::saturation_pressure_t (t);
            const auto
        d7 = r7::r2::density_pt (p, t);
            const auto
        [ d, info ] = r6_inverse::density_pt (p, t, d7, info::convergence);
        if (!info.converged) return 0.;
        return d;
    }
        auto
    d_sat_l_t (double t)
    {
        if (t < 273.16 || t > 647.096) return 0.;
            const auto
        p = r7::saturation_pressure_t (t);
            const auto
        [ d, info ] = r6_inverse::density_pt (p, t, 1000., info::convergence);
        if (!info.converged) return 0.;
        return d;
    }
        double
    d_sat_g;
        double
    d_sat_l;
public:
        void
    init (double t) override
    {
        d_sat_g = d_sat_g_t (t);
        d_sat_l = d_sat_l_t (t);
    };
        bool
    operator () (double d, double t) const override
    {
        return (t <= 273.16) || (t < 647.096 && d > d_sat_g && d < d_sat_l);
    };
};
    inline auto
exclusion_dt = exclusion_dt_t {};

    class
exclusion_tp_r7_t
    : public exclusion_base_t
{
public:
        bool
    operator () (double t, double p) const override
    {
        return t > 1074.15 && p > 50e6; 
    };
};
    inline auto
exclusion_tp_r7 = exclusion_tp_r7_t {};

    class
exclusion_tp_fluid_t
    : public exclusion_base_t
{
        double
    pmih;
        double
    pmiii;
        double
    pmv;
        double
    pmvi;
        double
    pmvii;
        double
    ps;
public:
        void
    init (double t) override
    {
        pmih  = r14::ih::melting_pressure_t  (t);
        pmiii = r14::iii::melting_pressure_t (t);
        pmv   = r14::v::melting_pressure_t   (t);
        pmvi  = r14::vi::melting_pressure_t  (t);
        pmvii = r14::vii::melting_pressure_t (t);
        ps    = r14::sublimation_pressure_t  (t);
    };
        bool
    operator () (double t, double p) const override
    {
        return
               (p <  611.657   && t <= 273.16 && p > ps)
            || (p >= 611.657   && p <  208.566e6 && t < 273.15  && p < pmih)
            || (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p > pmiii)
            || (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p > pmv)
            || (p >= 632.4e6   && p <  2216e6    && t < 355.    && p > pmvi)
            || (p >= 2216e6    && p <  2.06e10   && t < 715.    && p > pmvii)
            || (p >= 2.06e10   && t < 715.)
        ;
    };
};
    inline auto
exclusion_tp_fluid = exclusion_tp_fluid_t {};
    class
exclusion_tp_ice_t
    : public exclusion_tp_fluid_t
{
public:
        void
    init (double t) override
    {
        exclusion_tp_fluid_t::init (t);
    }
        bool
    operator () (double t, double p) const override
    {
        return !exclusion_tp_fluid_t::operator () (t, p);
    }
};
    inline auto
exclusion_tp_ice = exclusion_tp_ice_t {};
    class
exclusion_tp_r12_t
    : public exclusion_tp_fluid_t
{
public:
        void
    init (double t) override
    {
        exclusion_tp_fluid_t::init (t);
    }
        bool
    operator () (double t, double p) const override
    {
        return exclusion_tp_fluid_t::operator () (t, p) 
            || (p >  500e6 && t > 373.15)
            || (p >  350e6 && t > 433.15)
            || (p >  300e6 && t > 873.15)
            || t > 1173.15
        ;
    }
};
    inline auto
exclusion_tp_r12 = exclusion_tp_r12_t {};
    struct
range_t
{
        std::string
    name;
        std::pair <double, double>
    x;
        std::pair <double, double>
    y;
        std::pair <bool, bool>
    logscale = { false, false };
        exclusion_base_t&
    exclude;
        std::pair <int, int>
    values_count = { 501, 501 };
};

    class
exclusion_tp_supercooled_t
    : public exclusion_base_t
{
        double
    pmih;
        double
    pmiii;
        double
    pmv;
        double
    pmvi;
        double
    pmvii;
        double
    ps;
        double
    p_hin;
public:
        void
    init (double t) override
    {
        pmih  = r14::ih::melting_pressure_t  (t);
        pmiii = r14::iii::melting_pressure_t (t);
        pmv   = r14::v::melting_pressure_t   (t);
        pmvi  = r14::vi::melting_pressure_t  (t);
        pmvii = r14::vii::melting_pressure_t (t);
        ps    = r14::sublimation_pressure_t  (t);
        p_hin = g12::homogeneous_ice_nucleation_limit_temperature_low_t (t);
    };
        bool
    operator () (double t, double p) const override
    {
            const auto
        t_hin = g12::homogeneous_ice_nucleation_limit_temperature_high_p (p) ;
        return 
               (t < 181.4)
            || (t < 235.159 && p < 198.9e6 && p < p_hin)
            || (p > 198.9e6 && t < t_hin)
        /* TODO
            || (p < 611657  && t <= 273.16 && p > ps)
            || (p >= 611.657   && p <  208.566e6 && t < 273.15  && p > pmih)
            || (p >= 208.566e6 && p <  350.1e6   && t < 256.164 && p < pmiii)
            || (p >= 350.1e6   && p <  632.4e6   && t < 273.31  && p < pmv)
        */
        ;
    };
};
    inline auto
exclusion_tp_supercooled = exclusion_tp_supercooled_t {};

    template <class F>
    struct
graph_t
{
        std::string
    ztag;
        std::string
    xtag;
        std::string
    ytag;
        std::function <F>
    function;
        bool
    skip_zlog = false;
};

    template <class F>
    struct
topic_t
{
        std::string
    name;
        std::vector <range_t>
    ranges;
        std::vector <graph_t <F>>
    graphs;
        std::string
    title;
        std::string
    description;
        std::string
    substance;
};
    
    template <class F>
    auto
make_data_tp (topic_t <F> const& topic)
{
    for (auto&& range: topic.ranges)
    {
            auto
        o = std::unordered_map <std::string, std::ofstream> {};
        for (auto&& graph: topic.graphs)
        {
            o.emplace (
                  graph.ztag
                , graph.ztag + '_' + topic.name + '_' + range.name + ".dat"
            );
        }
            const auto
        dt = range.logscale.first 
            ? exp (log (range.x.second / range.x.first) / (range.values_count.first - 1))
            : (range.x.second - range.x.first) / (range.values_count.first - 1);
            const auto
        dp = range.logscale.second
            ? exp (log (range.y.second / range.y.first) / (range.values_count.second - 1))
            : (range.y.second - range.y.first) / (range.values_count.second - 1);
        for (auto it = 0; it < range.values_count.first; ++it)
        {
                const auto
            t = range.logscale.first 
                ? range.x.first * pow (dt, it)
                : range.x.first + dt * it;
            range.exclude.init (t);
            for (auto ip = 1; ip < range.values_count.second; ++ip)
            {
                    const auto
                p = range.logscale.second
                    ? range.y.first * pow (dp, ip)
                    : range.y.first + ip * dp;
                if (range.exclude (t, p)) continue;
                for (auto&& graph: topic.graphs)
                {
                    o.at (graph.ztag) 
                        << t 
                        << " " 
                        << p 
                        << " " << graph.function (t, p) << "\n"
                    ;
                }
            }
        }
    }
}
    template <class F>
    auto
make_data_dt (topic_t <F> const& topic)
{
    for (auto&& range: topic.ranges)
    {
            auto
        o = std::unordered_map <std::string, std::ofstream> {};
        for (auto&& graph: topic.graphs)
        {
            o.emplace (
                  graph.ztag
                , graph.ztag + '_' + topic.name + '_' + range.name + ".dat"
            );
        }
            const auto
        dd = range.logscale.first 
            ? exp (log (range.x.second / range.x.first) / (range.values_count.first - 1))
            : (range.x.second - range.x.first) / (range.values_count.first - 1);
            const auto
        dt = range.logscale.second
            ? exp (log (range.y.second / range.y.first) / (range.values_count.second - 1))
            : (range.y.second - range.y.first) / (range.values_count.second - 1);
        for (auto it = 0; it < range.values_count.second; ++it)
        {
                const auto
            t = range.logscale.second
                ? range.y.first * pow (dt, it)
                : range.y.first + dt * it
            ;
            range.exclude.init (t);
            for (auto id = 1; id < range.values_count.first; ++id)
            {
                    const auto
                d = range.logscale.first
                    ? range.x.first * pow (dd, id)
                    : range.x.first + id * dd
                ;
                if (range.exclude (d, t)) continue;
                for (auto&& graph: topic.graphs)
                {
                    o.at (graph.ztag) 
                        << d 
                        << " " 
                        << t 
                        << " " << graph.function (d, t) << "\n"
                    ;
                }
            }
        }
    }
}
// In the T-P plane, use r6_inverse to get D and use D-T functions.
    template <class F>
    auto
make_data_r6_inverse (topic_t <F> const& topic)
{
    for (auto&& range: topic.ranges)
    {
            auto
        o = std::unordered_map <std::string, std::ofstream> {};
        for (auto&& graph: topic.graphs)
        {
            o.emplace (
                  graph.ztag
                , graph.ztag + '_' + topic.name + '_' + range.name + ".dat"
            );
        }
            auto
        e = std::ofstream { topic.name + '_' + range.name + '_' + "convergence_failed.dat" };
            auto
        c = std::ofstream { topic.name + '_' + range.name + '_' + "_convergence_failed.txt" };
            auto
        j = std::ofstream { topic.name + '_' + range.name + '_' + "_convergence_iterations.txt" };
            const auto
        dt = range.logscale.first 
            ? exp (log (range.x.second / range.x.first) / (range.values_count.first - 1))
            : (range.x.second - range.x.first) / (range.values_count.first - 1);
            const auto
        dp = range.logscale.second
            ? exp (log (range.y.second / range.y.first) / (range.values_count.second - 1))
            : (range.y.second - range.y.first) / (range.values_count.second - 1);
        for (auto it = 0; it < range.values_count.first; ++it)
        {
                const auto
            t = range.logscale.first 
                ? range.x.first * pow (dt, it)
                : range.x.first + dt * it;
            range.exclude.init (t);
            for (auto ip = 1; ip < range.values_count.second; ++ip)
            {
                    const auto
                p = range.logscale.second
                    ? range.y.first * pow (dp, ip)
                    : range.y.first + ip * dp;
                if (range.exclude (t, p)) continue;
                    const auto
                di = (p > 100e6 || (t > 1073.15 && p > 50e6))
                    ? 1000.
                    : r7::density_pt (p, t)
                ;
                    const auto
                [ d, i ] = r6_inverse::density_pt (p, t, di, info::convergence);
                if (!i.converged)
                {
                    e << t << " " << p << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                for (auto&& graph: topic.graphs)
                {
                    o.at (graph.ztag) 
                        << t 
                        << " " 
                        << p 
                        << " " << graph.function (d, t) 
                        << "\n"
                    ;
                }
                j << t << " " << p << ' ' << i.convergence.size () << "\n";
                continue;
            }
        }
    }
}

    template <class F>
    auto
make_data_r7_vs_r6_inverse (topic_t <F> const& topic)
{
    for (auto&& range: topic.ranges)
    {
            auto
        o = std::unordered_map <std::string, std::ofstream> {};
        for (auto&& graph: topic.graphs)
        {
            o.emplace (
                  graph.ztag
                , graph.ztag + '_' + topic.name + '_' + range.name + ".dat"
            );
        }
            auto
        e = std::ofstream { topic.name + '_' + range.name + '_' + "convergence_failed.dat" };
            auto
        c = std::ofstream { topic.name + '_' + range.name + '_' + "_convergence_failed.txt" };
            auto
        j = std::ofstream { topic.name + '_' + range.name + '_' + "_convergence_iterations.txt" };
            const auto
        dt = range.logscale.first 
            ? exp (log (range.x.second / range.x.first) / (range.values_count.first - 1))
            : (range.x.second - range.x.first) / (range.values_count.first - 1);
            const auto
        dp = range.logscale.second
            ? exp (log (range.y.second / range.y.first) / (range.values_count.second - 1))
            : (range.y.second - range.y.first) / (range.values_count.second - 1);
        for (auto it = 0; it < range.values_count.first; ++it)
        {
                const auto
            t = range.logscale.first 
                ? range.x.first * pow (dt, it)
                : range.x.first + dt * it;
            range.exclude.init (t);
            for (auto ip = 1; ip < range.values_count.second; ++ip)
            {
                    const auto
                p = range.logscale.second
                    ? range.y.first * pow (dp, ip)
                    : range.y.first + ip * dp;
                if (range.exclude (t, p)) continue;
                    const auto
                d7 = r7::density_pt (p, t);
                    const auto
                [ d6, i ] = r6_inverse::density_pt (p, t, d7, info::convergence);
                if (!i.converged)
                {
                    e << t << " " << p << "\n";
                    c << format ("At {} K, {:e} Pa, with initial guess {} kg/m3\n", t, p, d7);
                    for (auto&& [ v, f ,d ]: i.convergence)
                    {
                        c << format ("  v: {}, f:{}, df: {}\n", v, f ,d);
                    }
                    continue;
                }
                for (auto&& graph: topic.graphs)
                {
                    o.at (graph.ztag) 
                        << t 
                        << " " 
                        << p 
                        << " " << graph.function (t, p, d6) 
                        << "\n"
                    ;
                }
                j << t << " " << p << ' ' << i.convergence.size () << "\n";
                continue;
            }
        }
    }
}

    template <class F>
    auto
make_gpl (topic_t <F> const& topic)
{
    for (auto&& range: topic.ranges)
    {
        for (auto&& graph: topic.graphs)
        {
                const auto
            base = graph.ztag + '_' + topic.name + '_' + range.name;
                auto
            o = std::ofstream (base + "_zlin.gpl");
            o << format (R"(set terminal pngcairo size 1280,960 enhanced
set xlabel "{}"
set ylabel "{}"
set cblabel "{}"
set ytics format "%.0e"
set cbtics format "%.0e"
set xrange [{}:{}]
set yrange [{}:{}]
set title "{}"
)"
                , quantities.at (graph.xtag).label + '[' + quantities.at (graph.xtag).unit + ']'
                , quantities.at (graph.ytag).label + '[' + quantities.at (graph.ytag).unit + ']'
                , quantities.at (graph.ztag).label + '[' + quantities.at (graph.ztag).unit + ']'
                , range.x.first  * quantities.at (graph.xtag).scale
                , range.x.second * quantities.at (graph.xtag).scale
                , range.y.first  * quantities.at (graph.ytag).scale
                , range.y.second * quantities.at (graph.ytag).scale
                , topic.substance + " " + quantities.at (graph.ztag).label + " (" + quantities.at (graph.xtag).label + ", " + quantities.at (graph.ytag).label + ')'
            );
            if (range.logscale.first) o << "set logscale x\n";
            if (range.logscale.second) o << "set logscale y\n";
            o << format (R"(set output "{}_zlin.png"
plot "{}.dat" using ($1*{}):($2*{}):($3*{}) w p pt 5 ps 0.3 lc palette z notitle
)"
                , base
                , base
                , quantities.at (graph.xtag).scale
                , quantities.at (graph.ytag).scale
                , quantities.at (graph.ztag).scale
            );
            o.close ();
            if (graph.skip_zlog) continue;
            o.open (base + "_zlog.gpl");
            o << format (R"(set terminal pngcairo size 1280,960 enhanced
set xlabel "{}"
set ylabel "{}"
set cblabel "{}"
set ytics format "%.0e"
set cbtics format "%.0e"
set xrange [{}:{}]
set yrange [{}:{}]
set title "{}"
set logscale cb
)"
                , quantities.at (graph.xtag).label + '[' + quantities.at (graph.xtag).unit + ']'
                , quantities.at (graph.ytag).label + '[' + quantities.at (graph.ytag).unit + ']'
                , quantities.at (graph.ztag).label + '[' + quantities.at (graph.ztag).unit + ']'
                , range.x.first  * quantities.at (graph.xtag).scale
                , range.x.second * quantities.at (graph.xtag).scale
                , range.y.first  * quantities.at (graph.ytag).scale
                , range.y.second * quantities.at (graph.ytag).scale
                , topic.substance + " " + quantities.at (graph.ztag).label + " (" + quantities.at (graph.xtag).label + ", " + quantities.at (graph.ytag).label + ')'
            );
            if (range.logscale.first) o << "set logscale x\n";
            if (range.logscale.second) o << "set logscale y\n";
            o << format (R"(set output "{}_zlog.png"
plot "{}.dat" using ($1*{}):($2*{}):($3*{}) w p pt 5 ps 0.3 lc palette z notitle
)"
                , base
                , base
                , quantities.at (graph.xtag).scale
                , quantities.at (graph.ytag).scale
                , quantities.at (graph.ztag).scale
            );
        }
    }
}
