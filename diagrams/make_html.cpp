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
        </p>
)";
        auto
    content = std::string {};
        auto
    toc = std::string { "<p>Contents:</p><ul>\n<li><a href='#lines'>Phases boundaries</a></li>\n" };
    for (auto&& topic: topics)
    {
        content += "<h2 id='" + topic.name + "'>" + topic.title + "</h2>\n<p>" + topic.description + "</p>\n";
        toc += "<li><a href='#" + topic.name + "'>" + topic.title + "</a></li>\n";
        for (auto&& graph: topic.graphs)
        {
            content += format (
                  "<h3>{} ({}, {})</h3>\n<div class='thumb-container'>\n"
                , quantities.at (graph.ztag).label
                , quantities.at (graph.xtag).label
                , quantities.at (graph.ytag).label
            );
            for (auto&& range: topic.ranges)
            {
                    const auto
                base = graph.ztag + '_' + topic.name + '_' + range.name;
                content += format (
                      "<div class='thumb'><a href='{}_zlin.png'><img src='{}_zlin-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "linear in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
                if (graph.skip_zlog) continue;
                content += format (
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
            content += format (
                  "<h3>{} ({}, {})</h3>\n<div class='thumb-container'>\n"
                , quantities.at (graph.ztag).label
                , quantities.at (graph.xtag).label
                , quantities.at (graph.ytag).label
            );
            for (auto&& range: topic.ranges)
            {
                    const auto
                base = graph.ztag + '_' + topic.name + '_' + range.name;
                content += format (
                      "<div class='thumb'><a href='{}_zlin.png'><img src='{}_zlin-thumb.png'/><div>"
                    , base
                    , base
                );
                if (range.logscale.first)  content += "logarithmic in " + quantities.at (graph.xtag).label + ", ";
                if (range.logscale.second) content += "logarithmic in " + quantities.at (graph.ytag).label + ", ";
                content += "linear in " + quantities.at (graph.ztag).label + "</div></a></div>\n";
                if (graph.skip_zlog) continue;
                content += format (
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
<div class='thumb-container'>
<div class='thumb'><a href='lines_tp.png'><img src='lines_tp-thumb.png'/><div>In the Temperature-Pressure plane</div></a></div>
<div class='thumb'><a href='lines_dt.png'><img src='lines_dt-thumb.png'/><div>In the Density-Temperature plane</div></a></div>
</div>
)";
    o 
        << content 
        << "<p>Updated " + format ("{:%Y-%m-%d}", std::chrono::system_clock::now ()) + "</p>"
        << "</body></html>"
    ;
}
