#include "diagram.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws;

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
                [ d6, i ] = r6_inverse::density_pt (p, t, d7, p * 1e-7, info::convergence);
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
    struct
exclusion_r7_vs_r6_inverse_t
    : exclusion_base_t
{
        void
    init (double) override
    {};
        bool
    operator () (double t, double p) const override
    {
        return t > 1074.15 && p > 50e6; 
    };
};
    auto
exclusion_r7_vs_r6_inverse = exclusion_r7_vs_r6_inverse_t {};

    double
yop (double d, double)
{
    return d;
}

    struct
delta
{
        std::function <double (double, double)>
    f6;
        std::function <double (double, double)>
    f7;
        template <class F6, class F7>
    delta (F6 a, F7 b)
        : f6 (a)
        , f7 (b)
    {}
        auto
    operator () (double t, double p, double d)
    {
            const auto
        a = f6 (d, t);
            const auto
        b = f7 (t, p);
        return fabs ((a - b) / a);
    }
};


    const auto
topic_r7_vs_r6_inverse = topic_t <double (double, double, double)>
{
      .name = "r7_vs_r6_inverse"
    , .ranges = {
          {
              "ylin"
            , { 273.16, 2273.16 }
            , { 0.,     100e6   }
            , { false,  false   }
            , exclusion_r7_vs_r6_inverse
          }
        , {
              "ylog"
            , { 273.16, 2273.16 }
            , { 1.    , 100e6   }
            , { false, true    }
            , exclusion_r7_vs_r6_inverse
          }
      }
    
    , .graphs = {
          { "delta_density",                              "temperature", "pressure", delta { yop,                                                          r7::density_tp                              <double, double> }}
        , { "delta_massic_enthalpy",                      "temperature", "pressure", delta { r6::massic_enthalpy_dt                      <double, double>, r7::massic_enthalpy_tp                      <double, double> }}
        , { "delta_massic_internal_energy",               "temperature", "pressure", delta { r6::massic_internal_energy_dt               <double, double>, r7::massic_internal_energy_tp               <double, double> }}
        , { "delta_massic_entropy",                       "temperature", "pressure", delta { r6::massic_entropy_dt                       <double, double>, r7::massic_entropy_tp                       <double, double> }}
        , { "delta_massic_isobaric_heat_capacity",        "temperature", "pressure", delta { r6::massic_isobaric_heat_capacity_dt        <double, double>, r7::massic_isobaric_heat_capacity_tp        <double, double> }}
        , { "delta_massic_isochoric_heat_capacity",       "temperature", "pressure", delta { r6::massic_isochoric_heat_capacity_dt       <double, double>, r7::massic_isochoric_heat_capacity_tp       <double, double> }}
        , { "delta_speed_of_sound",                       "temperature", "pressure", delta { r6::speed_of_sound_dt                       <double, double>, r7::speed_of_sound_tp                       <double, double> }}
        , { "delta_isobaric_cubic_expansion_coefficient", "temperature", "pressure", delta { r6::isobaric_cubic_expansion_coefficient_dt <double, double>, r7::isobaric_cubic_expansion_coefficient_tp <double, double> }}
        , { "delta_isothermal_compressibility",           "temperature", "pressure", delta { r6::isothermal_compressibility_dt           <double, double>, r7::isothermal_compressibility_tp           <double, double> }}
        , { "delta_relative_pressure_coefficient",        "temperature", "pressure", delta { r6::relative_pressure_coefficient_dt        <double, double>, r7::relative_pressure_coefficient_tp        <double, double> }}
        , { "delta_isothermal_stress_coefficient",        "temperature", "pressure", delta { r6::isothermal_stress_coefficient_dt        <double, double>, r7::isothermal_stress_coefficient_tp        <double, double> }}
      }
    , .title = "R7 vs. R6 inverse"
    , .description = R"(
<p>Comparison between the values obtained via R7 and R6 inverse.</p>
<p>Range: (173.15 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 100 MPa) and (1073.15 K &le; T &le; 2273.15 K and 0 MPa &le; P &le; 50 MPa)</p>
    )"
};


