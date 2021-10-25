#define ISTO_IAPWS_FORCE_RELAXED 1
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <fmt/format.h>
    using fmt::print, fmt::format;
    using namespace std::literals;
    struct
data_t
{
        std::string
    data_filename;
        std::string
    image_filename;
        std::string
    logscale;
        double
    xscale;
        double
    yscale;
        double
    zscale;
        std::array <double, 2>
    xrange;
        std::array <double, 2>
    yrange;
        std::string
    xlabel;
        std::string
    ylabel;
        std::string
    cblabel;
        std::string
    title;
        std::array <int, 2>
    size;
        std::array <int, 3>
    fields;
};

    struct
quantity_data_t
{
        std::string
    unit;
        double
    scale;
};


    auto const
quantities = std::unordered_map 
{
        std::pair
       { "Density"s,                   quantity_data_t                            { "kg/m^3"s,   1.,  }}
    ,  { "Pressure",                                                              { "MPa",       1e-6,}}
    ,  { "Temperature",                                                           { "K",         1.,  }}
    ,  { "Massic enthalpy",                                                       { "kJ/kg",     1e3, }}
    ,  { "Massic entropy",                                                        { "kJ/kg/K",   1e3, }}
    ,  { "Massic internal energy",                                                { "kJ/kg",     1e3, }}
    ,  { "Massic Gibbs free energy",                                              { "kJ/kg",     1e3, }}
    ,  { "Massic isobaric heat capacity",                                         { "kJ/kg/K",   1e3, }}
    ,  { "Massic isochoric heat capacity",                                        { "kJ/kg/K",   1e3, }}
    ,  { "Speed of sound",                                                        { "m/s",       1.,  }}
    ,  { "Isobaric cubic expansion coefficient",                                  { "K^{ -1}",   1.,  }}
    ,  { "Isothermal compressibility",                                            { "MPa^{ -1}", 1e6, }}
    ,  { "Isothermal stress coefficient",                                         { "kg/m^3",    1.,  }}
    ,  { "Relative pressure coefficient",                                         { "K^{ -1}",   1.,  }}
    ,  { "Relative pressure coefficient",                                         { "K^{ -1}",   1.,  }}
    ,  { "Iteration count",                                                       { "",          1.,  }}
    ,  { "Density difference, relative to R6 value",                              { "",          1.,  }}
    ,  { "Massic enthalpy difference, relative to R6 value",                      { "",          1.,  }}
    ,  { "Massic entropy difference, relative to R6 value",                       { "",          1.,  }}
    ,  { "Massic internal energy difference, relative to R6 value",               { "",          1.,  }}
    ,  { "Massic isobaric heat capacity difference, relative to R6 value",        { "",          1.,  }}
    ,  { "Massic isochoric heat capacity difference, relative to R6 value",       { "",          1.,  }}
    ,  { "Speed of sound difference, relative to R6 value",                       { "",          1.,  }}
    ,  { "Isobaric cubic expansion coefficient difference, relative to R6 value", { "",          1.,  }}
    ,  { "Isothermal compressibility difference, relative to R6 value",           { "",          1.,  }}
    ,  { "Isothermal stress coefficient difference, relative to R6 value",        { "",          1.,  }}
    ,  { "Relative pressure coefficient difference, relative to R6 value",        { "",          1.,  }}
};
    auto const
iapws_r7_graphs = std::vector
{
        std::tuple 
      { "density_r7"s,                             "Density"s,                             "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "massic_enthalpy_r7",                      "Massic enthalpy",                      "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_entropy_r7",                       "Massic entropy",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_internal_energy_r7",               "Massic internal energy",               "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isobaric_heat_capacity_r7",        "Massic isobaric heat capacity",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isochoric_heat_capacity_r7",       "Massic isochoric heat capacity",       "Temperature",  "Pressure",  1, 2, 3 }
    , { "speed_of_sound_r7",                       "Speed of sound",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "isobaric_cubic_expansion_coefficient_r7", "Isobaric cubic expansion coefficient", "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_compressibility_r7",           "Isothermal compressibility",           "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_stress_coefficient_r7",        "Isothermal stress coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "relative_pressure_coefficient_r7",        "Relative pressure coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
};
    auto const
iapws_r6_inv_graphs = std::vector
{
        std::tuple 
      { "density_r6_inv"s,                             "Density"s,                             "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "massic_enthalpy_r6_inv",                      "Massic enthalpy",                      "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_entropy_r6_inv",                       "Massic entropy",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_internal_energy_r6_inv",               "Massic internal energy",               "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_gibbs_free_energy_r6_inv",             "Massic Gibbs free energy",             "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isobaric_heat_capacity_r6_inv",        "Massic isobaric heat capacity",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isochoric_heat_capacity_r6_inv",       "Massic isochoric heat capacity",       "Temperature",  "Pressure",  1, 2, 3 }
    , { "speed_of_sound_r6_inv",                       "Speed of sound",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "isobaric_cubic_expansion_coefficient_r6_inv", "Isobaric cubic expansion coefficient", "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_compressibility_r6_inv",           "Isothermal compressibility",           "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_stress_coefficient_r6_inv",        "Isothermal stress coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "relative_pressure_coefficient_r6_inv",        "Relative pressure coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
};
    auto const
iapws_r6_inv_ext_graphs = std::vector
{
        std::tuple 
      { "density_r6_inv_ext"s,                             "Density"s,                             "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "massic_enthalpy_r6_inv_ext",                      "Massic enthalpy",                      "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_entropy_r6_inv_ext",                       "Massic entropy",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_internal_energy_r6_inv_ext",               "Massic internal energy",               "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_gibbs_free_energy_r6_inv_ext",             "Massic Gibbs free energy",             "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isobaric_heat_capacity_r6_inv_ext",        "Massic isobaric heat capacity",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isochoric_heat_capacity_r6_inv_ext",       "Massic isochoric heat capacity",       "Temperature",  "Pressure",  1, 2, 3 }
    , { "speed_of_sound_r6_inv_ext",                       "Speed of sound",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "isobaric_cubic_expansion_coefficient_r6_inv_ext", "Isobaric cubic expansion coefficient", "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_compressibility_r6_inv_ext",           "Isothermal compressibility",           "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_stress_coefficient_r6_inv_ext",        "Isothermal stress coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "relative_pressure_coefficient_r6_inv_ext",        "Relative pressure coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
};
    auto const
iapws_r7_vs_r6_inv_graphs = std::vector
{
        std::tuple 
      { "density_r7_vs_r6_inv_error"s,                             "Density difference, relative to R6 value"s,                             "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "density_r7_vs_r6_inv_iter"s,                              "Iteration count"s,                                                      "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "massic_enthalpy_r7_vs_r6_inv_error",                      "Massic enthalpy difference, relative to R6 value",                      "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_entropy_r7_vs_r6_inv_error",                       "Massic entropy difference, relative to R6 value",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_internal_energy_r7_vs_r6_inv_error",               "Massic internal energy difference, relative to R6 value",               "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isobaric_heat_capacity_r7_vs_r6_inv_error",        "Massic isobaric heat capacity difference, relative to R6 value",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isochoric_heat_capacity_r7_vs_r6_inv_error",       "Massic isochoric heat capacity difference, relative to R6 value",       "Temperature",  "Pressure",  1, 2, 3 }
    , { "speed_of_sound_r7_vs_r6_inv_error",                       "Speed of sound difference, relative to R6 value",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error", "Isobaric cubic expansion coefficient difference, relative to R6 value", "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_compressibility_r7_vs_r6_inv_error",           "Isothermal compressibility difference, relative to R6 value",           "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_stress_coefficient_r7_vs_r6_inv_error",        "Isothermal stress coefficient difference, relative to R6 value",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "relative_pressure_coefficient_r7_vs_r6_inv_error",        "Relative pressure coefficient difference, relative to R6 value",        "Temperature",  "Pressure",  1, 2, 3 }
};
/*
    auto const
iapws_r6_graphs = std::vector
{
        std::tuple 
      { "massic_enthalpy_r6"s,                     "Massic enthalpy"s,                     "Temperature"s, "Pressure"s, 1, 2, 3 }
    , { "massic_entropy_r6",                       "Massic entropy",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_internal_energy_r6",               "Massic internal energy",               "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_gibbs_free_energy_r6",             "Massic Gibbs free energy",             "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isobaric_heat_capacity_r6",        "Massic isobaric heat capacity",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "massic_isochoric_heat_capacity_r6",       "Massic isochoric heat capacity",       "Temperature",  "Pressure",  1, 2, 3 }
    , { "speed_of_sound_r6",                       "Speed of sound",                       "Temperature",  "Pressure",  1, 2, 3 }
    , { "isobaric_cubic_expansion_coefficient_r6", "Isobaric cubic expansion coefficient", "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_compressibility_r6",           "Isothermal compressibility",           "Temperature",  "Pressure",  1, 2, 3 }
    , { "isothermal_stress_coefficient_r6",        "Isothermal stress coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
    , { "relative_pressure_coefficient_r6",        "Relative pressure coefficient",        "Temperature",  "Pressure",  1, 2, 3 }
};
*/
    auto const
info_r6 = R"(
<h2>R6: The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use</h2>
<p>Provides direct calculations of the properties of ordinary water (i.e. water and stean) in terms of (density, temperature).</p>
<p></p>
<p>
Implements the following:
<ul>
<li>W. Wagner and A. Pruß , "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use", Journal of Physical and Chemical Reference Data 31, 387-535 (2002) <a href="https://doi.org/10.1063/1.1461829">https://doi.org/10.1063/1.1461829</a></li>
</ul>
</p>
)";

    auto const
info_r6_inv = R"(
<h2>R6 inverse</h2>
<p>Provides indirect calculations of the properties of ordinary water in terms 
of (pressure, temperature) by inverting the P(&rho;, T) relation of R6.</p>
<p>Range: 173.15 K &le; T &le; 1273.15 K and 0 MPa &le; P &le; 1000 MPa </p>
)";

    auto const
info_r6_inv_ext = R"(
<h2>R6 inverse, extended to 5000 K, 100 GPa</h2>
<p>This is R6 inverse over an extended (P, T) range.</p>
<p>Range: 173.15 K &le; T &le; 5073.15 K and 0 MPa &le; P &le; 100 GPa </p>
)";

    auto const
info_r7 = R"(
<h2>R7: The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam</h2>
<p>Provides direct calculations of the properties of ordinary water (i.e. water and stean) in terms of (pressure, temperature).</p>
<p>It is an approximation of R6.</p>
<p>Range: (173.15 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 100 MPa) and (1073.15 K &le; T &le; 2273.15 K and 0 MPa &le; P &le; 50 MPa)</p>
<p>
Implements the following:
<ul>
<li>Wagner, W., Cooper, J. R., Dittmann, A., Kijima, J., Kretzschmar, H., Kruse, A., Mareš, R., Oguchi, K., Sato, H., Stöcker, I., Sǐfner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., and Willkommen, T. (January 1, 2000). "The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam ." ASME. J. Eng. Gas Turbines Power. January 2000; 122(1): 150–184. <a href="https://doi.org/10.1115/1.483186">https://doi.org/10.1115/1.483186</a></li>
<li>Kretzschmar, H., Cooper, J. R., Dittmann, A., Friend, D. G., Gallagher, J. S., Knobloch, K., Mareš, R., Miyagawa, K., Stöcker, I., Trübenbach, J., Wagner, W., and Willkommen, T. (June 22, 2004). "Supplementary Backward Equations for Pressure as a Function of Enthalpy and Entropy li(h,s) to the Industrial Formulation IAPWS-IF97 for Water and Steam." ASME. J. Eng. Gas Turbines Power. July 2006; 128(3): 702–713. <a href="https://doi.org/10.1115/1.1915392">https://doi.org/10.1115/1.1915392</a></li>
<li>Kretzschmar, H., Cooper, J. R., Dittmann, A., Friend, D. G., Gallagher, J. S., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Stöcker, I., Wagner, W., and Weber, I. (January 10, 2006). "Supplementary Backward Equations T(li,h)⁠, v(li,h)⁠, and T(li,s)⁠, v(li,s) for the Critical and Supercritical Regions (Region 3) of the Industrial Formulation IAPWS-IF97 for Water and Steam." ASME. J. Eng. Gas Turbines Power. January 2007; 129(1): 294–303. <a href="https://doi.org/10.1115/1.2181598">https://doi.org/10.1115/1.2181598</a></li>
<li>Kretzschmar, H., Cooper, J. R., Gallagher, J. S., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Span, R., Stöcker, I., Wagner, W., and Weber, I. (January 16, 2007). "Supplementary Backward Equations li(h,s) for the Critical and Supercritical Regions (Region 3), and Equations for the Two-Phase Region and Region Boundaries of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam." ASME. J. Eng. Gas Turbines Power. October 2007; 129(4): 1125–1137. <a href="https://doi.org/10.1115/1.2719267">https://doi.org/10.1115/1.2719267</a></li>
<li>Kretzschmar, H., Harvey, A. H., Knobloch, K., Mareš, R., Miyagawa, K., Okita, N., Span, R., Stöcker, I., Wagner, W., and Weber, I. (April 13, 2009). "Supplementary Backward Equations v(li,T) for the Critical and Supercritical Regions (Region 3) of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam." ASME. J. Eng. Gas Turbines Power. July 2009; 131(4): 043101. <a href="https://doi.org/10.1115/1.3028630">https://doi.org/10.1115/1.3028630</a></li>
</ul>
</p>
)";

    auto const
info_r7_vs_r6 = R"(
<h2>R7 vs. R6 inverse</h2>
<p>Comparison between the values obtained via R7 and R6 inverse.</p>
<p>Range: (173.15 K &le; T &lt; 1273.15 K and 0 MPa &le; P &le; 100 MPa) and (1073.15 K &le; T &le; 2273.15 K and 0 MPa &le; P &le; 50 MPa)</p>
)";

    auto const
model_list = 
{
      "iapws_r7"
    , "iapws_r6_inverse"
    , "iapws_r6_inverse_extended"
    , "iapws_r7_vs_r6_inverse"
};

    auto const
model_name = std::unordered_map
{
        std::pair
        { "iapws_r7",                   "R7" }
    ,   { "iapws_r6_inverse",           "R6 inverse" }
    ,   { "iapws_r6_inverse_extended",  "R6 inverse, extended range" }
    ,   { "iapws_r7_vs_r6_inverse",     "Comparison between R7 and R6 inverse" }
};

    auto const
models = std::unordered_map
{
        std::pair
      { "iapws_r7",       std::tuple { std::array { 273., 2273. }, std::array { 0., 1e2 }, std::array { 1e-6, 1e2 }, iapws_r7_graphs,           info_r7         }}
    , { "iapws_r6_inverse",          { std::array { 173., 1273. }, std::array { 0., 1e3 }, std::array { 1e-6, 1e3 }, iapws_r6_inv_graphs,       info_r6_inv     }}
    , { "iapws_r6_inverse_extended", { std::array { 173., 5073. }, std::array { 0., 1e5 }, std::array { 1e-6, 1e5 }, iapws_r6_inv_ext_graphs,   info_r6_inv_ext }}
    , { "iapws_r7_vs_r6_inverse",    { std::array { 273., 2273. }, std::array { 0., 1e2 }, std::array { 1e-6, 1e2 }, iapws_r7_vs_r6_inv_graphs, info_r7_vs_r6   }}
    //, { "IAPWS R6"s,                   { std::array { 173., 1273. }, std::array { 0., 1e3 }, std::array { 1e-6, 1e3 }, iapws_r6_graphs,                         }}
};

    auto
gen_data ()
{
        auto
    d = std::vector <data_t> {};
    for (auto&& [ model, A ]: models)
    {
            auto&&
        [ xrange, yrange1, yrange2, graphs, info ] = A;
        for (auto&& [ datafile, zi, xi, yi, f1, f2, f3 ]: graphs)
        {
                auto&&
            xd = quantities.at (xi);
                auto&&
            yd = quantities.at (yi);
                auto&&
            zd = quantities.at (zi);
            d.push_back ({
                  datafile + "_ylin"
                , datafile + "_ylin_zlin"
                , "" 
                , xd.scale
                , yd.scale
                , zd.scale
                , xrange
                , yrange1
                , xi + " [" + xd.unit + "]"
                , yi + " [" + yd.unit + "]"
                , zi + " [" + zd.unit + "]"
                , zi + " (" + xi + ", " + yi + ") according to " + model
                , { 1280 , 960 }
                , { f1, f2, f3 }
            });
            d.push_back ({
                  datafile + "_ylin"
                , datafile + "_ylin_zlog"
                , "cb" 
                , xd.scale
                , yd.scale
                , zd.scale
                , xrange
                , yrange1
                , xi + " [" + xd.unit + "]"
                , yi + " [" + yd.unit + "]"
                , zi + " [" + zd.unit + "]"
                , zi + " (" + xi + ", " + yi + ") according to " + model
                , { 1280 , 960 }
                , { f1, f2, f3 }
            });
            d.push_back ({
                  datafile + "_ylog"
                , datafile + "_ylog_zlin"
                , "y" 
                , xd.scale
                , yd.scale
                , zd.scale
                , xrange
                , yrange2
                , xi + " [" + xd.unit + "]"
                , yi + " [" + yd.unit + "]"
                , zi + " [" + zd.unit + "]"
                , zi + " (" + xi + ", " + yi + ") according to " + model
                , { 1280 , 960 }
                , { f1, f2, f3 }
            });
            d.push_back ({
                  datafile + "_ylog"
                , datafile + "_ylog_zlog"
                , "ycb" 
                , xd.scale
                , yd.scale
                , zd.scale
                , xrange
                , yrange2
                , xi + " [" + xd.unit + "]"
                , yi + " [" + yd.unit + "]"
                , zi + " [" + zd.unit + "]"
                , zi + " (" + xi + ", " + yi + ") according to " + model
                , { 1280 , 960 }
                , { f1, f2, f3 }
            });
        }
    }
    return d;
}
    auto
gen_gpl ()
{
        auto const
    data = gen_data ();
    for (auto&& d: data)
    {
            auto
        f = std::ofstream { d.image_filename + ".gpl" };
        f << format (R"(reset
set terminal pngcairo size {0},{1}
set xlabel "{2}"
set ylabel "{3}"
set cblabel "{4}"
set ytics format "%.0e"
set cbtics format "%.0e"
set xrange [{5}:{6}]
set yrange [{7}:{8}]
set title "{9}"
{10}
set output "{11}.png"
plot "{12}.dat" using (${16}*{13}):(${17}*{14}):(${18}*{15}) w p pt 5 ps 0.3 lc palette z notitle \
    , "sublimation.dat"  using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "white" notitle \
    , "saturation.dat"   using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "white" notitle \
    , "melting_ih.dat"   using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "black" notitle \
    , "melting_iii.dat"  using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "black" notitle \
    , "melting_v.dat"    using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "black" notitle \
    , "melting_vi.dat"   using ($1*{13}):($2*{14}) w l lw 3 lc rgbcolor "black" notitle

# vim:ft=gnuplot
    )"
            , d.size[0]
            , d.size[1]
            , d.xlabel
            , d.ylabel
            , d.cblabel
            , d.xrange[0]
            , d.xrange[1]
            , d.yrange[0]
            , d.yrange[1]
            , d.title
            , d.logscale.empty () ? "unset logscale"s : "set logscale "s + d.logscale
            , d.image_filename
            , d.data_filename
            , d.xscale
            , d.yscale
            , d.zscale
            , d.fields.at (0)
            , d.fields.at (1)
            , d.fields.at (2)
        );
    }
}
    auto
gen_html ()
{
        auto
    f = std::ofstream { "index.html" };
    f << R"( <!doctype html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <style type="text/css">
            body {
                font-family: sans-serif;
                color: rgb(102, 102, 102);
                max-width: 1360px;
                margin: auto;
            }
            a {
                color: rgb(51, 122, 183);
                text-decoration: none;
            }
            div.paragraph {
                margin-bottom: 60px;
            }
            nav {
                background-color: #d2722b;
                border-color: #e7e7e7;
                top: 0;
                border-width: 0 0 1px;
                min-height: 40px;
                margin-bottom: 40px;
                display: block;
            }
            nav div {
            }
            nav div a {
                color: white;
                background-color: transparent;
                text-decoration: none;
                float: left;
                height: 40px;
                margin: 10px 20px;
                font-size: 18px;
                line-height: 20px;
            }
            li {
                margin: 5px;
            }
            .thumb-container {
                display: grid;
                grid-template-columns: 1fr 1fr 1fr 1fr
            }
            .thumb {
                width: 320px;
                display: inline flow-root;
                margin: 10px;
            }
        </style>
        <link rel="icon" href="/cropped-favicon-32x32.jpg" sizes="32x32">
        <title>Water phase diagrams</title>
    </head>
    <body> 
        <p>/!\ Work in progress!</p>
        <h1>Water phase diagrams, according to IAPWS</h1>

        <h2>Phase boundaries</h2>
        <div class="thumb-container">
        <div class="thumb"><a href="lines.png"><img src="lines-thumb.png"/><div>Phase boundaries in the T-P plane</div></a></div>
        <div class="thumb"><a href="lines_dt.png"><img src="lines_dt-thumb.png"/><div>Phase boundaries in the &rho;-T plane</div></a></div>
        </div>)";
    //for (auto&& [ model, A ]: models)
    f << "<p>Content:</p><ul>\n";
    for (auto&& id: model_list)
    {
            auto&&
        n = model_name.at (id);
        f << format (R"(<li><a href="#{}">{}</a></li>)", id, n);
    }
    f << "</ul>\n";
    for (auto&& m: model_list)
    {
            auto&&
        [ xrange, yrange1, yrange2, graphs, info ] = models.at (m);
        f << format (R"(<span id="{}"/>{})", m, info);
        for (auto&& [ file, zi, xi, yi, f1, f2, f3 ]: graphs)
        {
            f << format (R"(            
            <h3>{0} ({1}, {2})</h3>
            <div class="thumb-container">
                <div class="thumb"><a href="{3}_ylin_zlin.png"><img src="{3}_ylin_zlin-thumb.png"/><div>Linear in pressure, linear in {0}</div></a></div>
                <div class="thumb"><a href="{3}_ylin_zlog.png"><img src="{3}_ylin_zlog-thumb.png"/><div>Linear in pressure, logarithmic in {0}</div></a></div>
                <div class="thumb"><a href="{3}_ylog_zlin.png"><img src="{3}_ylog_zlin-thumb.png"/><div>logarithmic in pressure, linear in {0}</div></a></div>
                <div class="thumb"><a href="{3}_ylog_zlog.png"><img src="{3}_ylog_zlog-thumb.png"/><div>logarithmic in pressure, logarithmic in {0}</div></a></div>
            </div>
)"
                , zi
                , xi
                , yi
                , file
                , info
            );
        }
    }
    f << "\n    </body>\n</html>";
}
    int
main ()
{
    gen_gpl ();
    {
            auto
        f = std::ofstream { "lines.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960
set xlabel "temperature [K]"
set ylabel "pressure [MPa]"
set cblabel "density [kg/m3]"
set ytics format "%.0e"
set cbtics format "%.0e"
set logscale y
set output "lines.png"
set key bottom
set title "Phases boundaries in the T-P plane"
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
      "sublimation.dat" using 1:($2*1e-6) w l t "sublimation line" \
    , "saturation.dat"  using 1:($2*1e-6) w l t "saturation line"  \
    , "melting_ih.dat"  using 1:($2*1e-6) w l t "melting line Ih"  \
    , "melting_iii.dat" using 1:($2*1e-6) w l t "melting line III" \
    , "melting_v.dat"   using 1:($2*1e-6) w l t "melting line V"   \
    , "melting_vi.dat"  using 1:($2*1e-6) w l t "melting line VI"  \
    , "melting_vii.dat" using 1:($2*1e-6) w l t "melting line VII"

# vim:ft=gnuplot
    )";
    }
    {
            auto
        f = std::ofstream { "lines_dt.gpl" };
        f << R"(reset
set terminal pngcairo size 1280,960
set ylabel "temperature [K]"
set xlabel "density [kg/m3]"
set output "lines_dt.png"
set key bottom
set title "Phases boundaries in the \U+03C1-T plane"
set xrange [-100:1900]
set yrange [200:700]
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
    , "melting_vii_dt.dat"    w l t "Melting line VII"

# vim:ft=gnuplot
    )";
    }
    gen_html ();
    return 0;
}

