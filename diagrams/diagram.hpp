#pragma once
#define ISTO_IAPWS_FORCE_RELAXED 1
#include <functional>
#include <string>
#include <unordered_map>
#include <fstream>
#include <fmt/format.h>
    using fmt::print, fmt::format;
    using namespace std::literals;

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
    ,  { "pressure",                                   { "Pressure",                                                              "MPa",       1e-6,}}
    ,  { "temperature",                                { "Temperature",                                                           "K",         1.,  }}
    ,  { "massic_enthalpy",                            { "Massic enthalpy",                                                       "kJ/kg",     1e3, }}
    ,  { "massic_entropy",                             { "Massic entropy",                                                        "kJ/kg/K",   1e3, }}
    ,  { "massic_internal_energy",                     { "Massic internal energy",                                                "kJ/kg",     1e3, }}
    ,  { "massic_gibbs_free_energy",                   { "Massic Gibbs free energy",                                              "kJ/kg",     1e3, }}
    ,  { "massic_isobaric_heat_capacity",              { "Massic isobaric heat capacity",                                         "kJ/kg/K",   1e3, }}
    ,  { "massic_isochoric_heat_capacity",             { "Massic isochoric heat capacity",                                        "kJ/kg/K",   1e3, }}
    ,  { "speed_of_sound",                             { "Speed of sound",                                                        "m/s",       1.,  }}
    ,  { "isobaric_cubic_expansion_coefficient",       { "Isobaric cubic expansion coefficient",                                  "K^{ -1}",   1.,  }}
    ,  { "isothermal_compressibility",                 { "Isothermal compressibility",                                            "MPa^{ -1}", 1e6, }}
    ,  { "isothermal_stress_coefficient",              { "Isothermal stress coefficient",                                         "kg/m^3",    1.,  }}
    ,  { "relative_pressure_coefficient",              { "Relative pressure coefficient",                                         "K^{ -1}",   1.,  }}
    ,  { "relative_pressure_coefficient",              { "Relative pressure coefficient",                                         "K^{ -1}",   1.,  }}
    ,  { "iteration_count",                            { "Iteration count",                                                       "",          1.,  }}
    ,  { "delta_density",                              { "Density difference, relative to R6 value",                              "",          1.,  }}
    ,  { "delta_massic_enthalpy",                      { "Massic enthalpy difference, relative to R6 value",                      "",          1.,  }}
    ,  { "delta_massic_entropy",                       { "Massic entropy difference, relative to R6 value",                       "",          1.,  }}
    ,  { "delta_massic_internal_energy",               { "Massic internal energy difference, relative to R6 value",               "",          1.,  }}
    ,  { "delta_massic_isobaric_heat_capacity",        { "Massic isobaric heat capacity difference, relative to R6 value",        "",          1.,  }}
    ,  { "delta_massic_isochoric_heat_capacity",       { "Massic isochoric heat capacity difference, relative to R6 value",       "",          1.,  }}
    ,  { "delta_speed_of_sound",                       { "Speed of sound difference, relative to R6 value",                       "",          1.,  }}
    ,  { "delta_isobaric_cubic_expansion_coefficient", { "Isobaric cubic expansion coefficient difference, relative to R6 value", "",          1.,  }}
    ,  { "delta_isothermal_compressibility",           { "Isothermal compressibility difference, relative to R6 value",           "",          1.,  }}
    ,  { "delta_isothermal_stress_coefficient",        { "Isothermal stress coefficient difference, relative to R6 value",        "",          1.,  }}
    ,  { "delta_relative_pressure_coefficient",        { "Relative pressure coefficient difference, relative to R6 value",        "",          1.,  }}
};

    struct
exclusion_base_t
{
        virtual void
    init (double)
    {}
        virtual bool
    operator () (double, double) const
    {
        return false;
    }
};

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
        std::function <F>//std::function <double (double, double)>
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
                , quantities.at (graph.ztag).label + " (" + quantities.at (graph.xtag).label + ", " + quantities.at (graph.ytag).label + ')'
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
                , quantities.at (graph.ztag).label + " (" + quantities.at (graph.xtag).label + ", " + quantities.at (graph.ytag).label + ')'
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
