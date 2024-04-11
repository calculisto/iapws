#include "diagram_all.hpp"
#include <chrono>
#include <fmt/chrono.h>

    int 
main ()
{
        auto
    o = std::ofstream { "index.html" };
    o << R"(<!doctype html>
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
        <p>/!\ Work in progress, you might want to press Ctrl-F5 from time to time</p>
        <h1>Phase diagrams of ordinary water</h1>
        <p>
        On this page are presented some phase diagrams for water. Every data point 
        plotted on the graphs shown below was computed using one or a combination of 
        several models published by the International Association for the 
        Properties of Water and Steam (IAPWS). These publications (Realeases in 
        the IAPWS lingo) are identified by the letter R and a number. Very 
        briefly:
        <ul>
        <li>R6 is the "master equation" that has been carefully fitted over a 
        wide range of experimental results and represents the state of the art 
        for the evaluation of the thermodynamic properties of water. 
        It is expressed in terms of Helmholtz energy as a function of density 
        and temperature, from which all the other properties can be obtained by 
        some combinations of derivatives.</li>
        <li>R7 is an approximation of R6 that allow direct and fast calculations 
        of properties, e.g. as functions of temperature and pressure, where R6 
        would require a computational expensive inversion. Below are presented 
        diagrams built with R7 as well as the same obtained by inverting R6 
        (under the sobriquet "R6 inverse"), and the differences between the two
        ("R7 vs. R6 inverse")</li>
        <li>R14 describes the boundaries between fluid water and various phases 
        of ice (a.k.a. melting and sublimation curves). The saturation curve can
        be directly calculated in the P-T plane using R7.</li>
        <li>R10 describes ice Ih.</li>
        <li>R12 provides the viscosity of liquid and gaseous water.</li>
        </ul>

        The implementation of the IAPWS models used to compute the data is 
        <a href="https://github.com/le-migou/iapws">isto/iapws</a>.

        When needed, the inversions were performed using the 
        <a href="https://github.com/le-migou/root_finding">isto/root_finding</a> 
        library.

        The plots were rendered using <a href="http://www.gnuplot.info/">Gnuplot</a>.
        </p>
        <p>
        Contact:
        <span style="unicode-bidi:bidi-override; direction: rtl;">
        rf.srnc@gnort-el.leunamme
        </span></span>
        </p>
)";
        auto
    content = std::string {};
        auto
    toc = std::string { "<p>Contents:</p><ul>\n<li><a href='#lines'>Phases boundaries</a></li>\n<li><a href='#spinodal'>Spinodal lines</a></li>\n" };
    for (auto&& topic: topics)
    {
        content += "<h2 id='" + topic.name + "'>" + topic.title + "</h2>\n<p>" + topic.description + "</p>\n";
        toc += "<li><a href='#" + topic.name + "'>" + topic.title + "</a></li>\n";
        for (auto&& graph: topic.graphs)
        {
            content += fmt::format (
                  "<h3>{} ({}, {})</h3>\n<div class='thumb-container'>\n"
                , quantities.at (graph.ztag).label
                , quantities.at (graph.xtag).label
                , quantities.at (graph.ytag).label
            );
            for (auto&& range: topic.ranges)
            {
                    const auto
                base = graph.ztag + '_' + topic.name + '_' + range.name;
                content += fmt::format (
                      "<div class='thumb'><a href='{}_zlin.png'><img src='{}_zlin-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "linear in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
                if (graph.skip_zlog) continue;
                content += fmt::format (
                      "<div class='thumb'><a href='{}_zlog.png'><img src='{}_zlog-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "logarithmic in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
            }
            content += "</div>\n";
        }
    }
    {
            auto const&
        topic = topic_r7_vs_r6_inverse;
        content += "<h2 id='" + topic.name + "'>" + topic.title + "</h2>\n<p>" + topic.description + "</p>\n";
        toc += "<li><a href='#" + topic.name + "'>" + topic.title + "</a></li>\n";
        for (auto&& graph: topic.graphs)
        {
            content += fmt::format (
                  "<h3>{} ({}, {})</h3>\n<div class='thumb-container'>\n"
                , quantities.at (graph.ztag).label
                , quantities.at (graph.xtag).label
                , quantities.at (graph.ytag).label
            );
            for (auto&& range: topic.ranges)
            {
                    const auto
                base = graph.ztag + '_' + topic.name + '_' + range.name;
                content += fmt::format (
                      "<div class='thumb'><a href='{}_zlin.png'><img src='{}_zlin-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "linear in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
                if (graph.skip_zlog) continue;
                content += fmt::format (
                      "<div class='thumb'><a href='{}_zlog.png'><img src='{}_zlog-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "logarithmic in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
            }
            content += "</div>\n";
        }
    }
    toc += "</ul>\n";
    o << toc;
    o << R"(
<h2 id='lines'>Phases boundaries</h2>
<p>Phases boundaries calculated using
<ul>
<li>Wagner, W., Cooper, J. R., Dittmann, A., Kijima, J., Kretzschmar, H., Kruse, A., Mareš, R., Oguchi, K., Sato, H., Stöcker, I., Sǐfner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., and Willkommen, T. (January 1, 2000). "The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam ." ASME. J. Eng. Gas Turbines Power. January 2000; 122(1): 150–184. <a href="https://doi.org/10.1115/1.483186">https://doi.org/10.1115/1.483186</a></li>
<li>Wolfgang Wagner, Thomas Riethmann, Rainer Feistel, Allan H. Harvey, "New Equations for the Sublimation Pressure and Melting Pressure of H2O Ice Ih", Journal of Physical and Chemical Reference Data 40, 043103 (2011) <a href="https://doi.org/10.1063">https://doi.org/10.1063</a></li>
</ul>
</p>
<div class='thumb-container'>
<div class='thumb'><a href='lines_tp.png'><img src='lines_tp-thumb.png'/><div>In the Temperature-Pressure plane</div></a></div>
<div class='thumb'><a href='lines_dt.png'><img src='lines_dt-thumb.png'/><div>In the Density-Temperature plane</div></a></div>
</div>
<h2 id='spinodal'>Spinodal lines</h2>
<p>
Here we use IAPWS 95 (a.k.a. R6) to find the spinodal lines, i.e. locus of points where the compressibility is zero. 
"Spinodal line (liq&#41;" uses R6 on the liquid side (densities &gt; critical density).
"Spinodal line (gas), candidate 1" and "Spinodal line (gas), candidate 1" use R6 on the vapor side (densities &lt; critical density).
It so happens that there are two solutions on this side, although the second one has probably no thermophysical sense. 
"Spinodal line (gas), R6 gas" uses the gas equation of state found in the paper of W. Wagner and A. Pruß (eq.&nbsp;3.2)
(but not on the official IAPWS guideline). Unfortunately, its validity range does not extend beyond a density of 55&nbsp;kg/m<spu>3</sup>


</p>
<div class='thumb-container'>
<div class='thumb'><a href='spinodal_tp.png'><img src='spinodal_tp-thumb.png'/><div>In the Temperature-Pressure plane</div></a></div>
<div class='thumb'><a href='spinodal_dt.png'><img src='spinodal_dt-thumb.png'/><div>In the Density-Temperature plane</div></a></div>
</div>
)";
    o 
        << content 
        << "<p>Updated " + fmt::format ("{:%Y-%m-%d}", std::chrono::system_clock::now ()) + "</p>"
        << "</body></html>"
    ;
}
