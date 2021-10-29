#pragma once
#include "diagram.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
#include "../include/isto/iapws/r14.hpp"
    using namespace isto::iapws;

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
                [ d, i ] = r6_inverse::density_pt (p, t, di, p * 1e-7, info::convergence);
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


    struct
exclusion_r6_inverse_t
    : exclusion_base_t
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
    auto
exclusion_r6_inverse = exclusion_r6_inverse_t {};
    const auto
topic_r6_inverse = topic_t <double (double, double)>
{
      .name = "r6_inverse"
    , .ranges = {
          {
              "ylin"
            , { 173.16, 1273.16 }
            , { 0.,     1000e6  }
            , { false,  false   }
            , exclusion_r6_inverse
          }
        , {
              "ylog"
            , { 173.16, 2273.16 }
            , { 1.    , 1000e6  }
            , { false,  true    }
            , exclusion_r6_inverse
          }
      }
    
    , .graphs = {
          { "density",                              "temperature", "pressure", [](double d, double) { return d; }                           }
        , { "massic_enthalpy",                      "temperature", "pressure", r6::massic_enthalpy_dt                      <double, double> }
        , { "massic_internal_energy",               "temperature", "pressure", r6::massic_internal_energy_dt               <double, double> }
        , { "massic_gibbs_free_energy",             "temperature", "pressure", r6::massic_gibbs_free_energy_dt             <double, double> }
        , { "massic_entropy",                       "temperature", "pressure", r6::massic_entropy_dt                       <double, double> }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", r6::massic_isobaric_heat_capacity_dt        <double, double> }
        , { "massic_isochoric_heat_capacity",       "temperature", "pressure", r6::massic_isochoric_heat_capacity_dt       <double, double> }
        , { "speed_of_sound",                       "temperature", "pressure", r6::speed_of_sound_dt                       <double, double> }
        , { "isobaric_cubic_expansion_coefficient", "temperature", "pressure", r6::isobaric_cubic_expansion_coefficient_dt <double, double> }
        , { "isothermal_compressibility",           "temperature", "pressure", r6::isothermal_compressibility_dt           <double, double> }
        , { "relative_pressure_coefficient",        "temperature", "pressure", r6::relative_pressure_coefficient_dt        <double, double> }
        , { "isothermal_stress_coefficient",        "temperature", "pressure", r6::isothermal_stress_coefficient_dt        <double, double> }
      }
    , .title = "R6 inverse"
    , .description = R"(
<p>Provides indirect calculations of the properties of ordinary water in terms 
of (pressure, temperature) by inverting the P(&rho;, T) relation of R6.</p>
<p>Range: 173.15 K &le; T &le; 1273.15 K and 0 MPa &le; P &le; 1000 MPa </p>
    )"
};

    auto
exclusion_r6_inverse_extended = exclusion_r6_inverse_t {};
    const auto
topic_r6_inverse_extended = topic_t <double (double, double)>
{
      .name = "r6_inverse_extended"
    , .ranges = {
          {
              "ylin"
            , { 173.16, 5073.16 }
            , { 0.,     100e9   }
            , { false,  false   }
            , exclusion_r6_inverse_extended
          }
        , {
              "ylog"
            , { 173.16, 5073.16 }
            , { 1.    , 100e9   }
            , { false,  true    }
            , exclusion_r6_inverse_extended
          }
      }
    
    , .graphs = {
          { "density",                              "temperature", "pressure", [](double d, double) { return d; }                           }
        , { "massic_enthalpy",                      "temperature", "pressure", r6::massic_enthalpy_dt                      <double, double> }
        , { "massic_internal_energy",               "temperature", "pressure", r6::massic_internal_energy_dt               <double, double> }
        , { "massic_gibbs_free_energy",             "temperature", "pressure", r6::massic_gibbs_free_energy_dt             <double, double> }
        , { "massic_entropy",                       "temperature", "pressure", r6::massic_entropy_dt                       <double, double> }
        , { "massic_isobaric_heat_capacity",        "temperature", "pressure", r6::massic_isobaric_heat_capacity_dt        <double, double> }
        , { "massic_isochoric_heat_capacity",       "temperature", "pressure", r6::massic_isochoric_heat_capacity_dt       <double, double> }
        , { "speed_of_sound",                       "temperature", "pressure", r6::speed_of_sound_dt                       <double, double> }
        , { "isobaric_cubic_expansion_coefficient", "temperature", "pressure", r6::isobaric_cubic_expansion_coefficient_dt <double, double> }
        , { "isothermal_compressibility",           "temperature", "pressure", r6::isothermal_compressibility_dt           <double, double> }
        , { "relative_pressure_coefficient",        "temperature", "pressure", r6::relative_pressure_coefficient_dt        <double, double> }
        , { "isothermal_stress_coefficient",        "temperature", "pressure", r6::isothermal_stress_coefficient_dt        <double, double> }
      }
    , .title = "R6 inverse, extended to 5000 K, 100 GPa"
    , .description = R"(
<p>This is R6 inverse over an extended (P, T) range.</p>
<p>The phase boundary for P &gt; 20 GPa is totally arbitrary, as no model is available in this range (AFAIK).</p>
<p>Range: 173.15 K &le; T &le; 5073.15 K and 0 MPa &le; P &le; 100 GPa </p>
      )"
};
