#include "diagram_r6.hpp"
#include "diagram_r6_inverse.hpp"
#include "diagram_r7.hpp"
#include "diagram_r7_vs_r6_inverse.hpp"

    const auto
topics = std::vector
{
      topic_r6
    , topic_r6_inverse
    , topic_r6_inverse_extended
    , topic_r7
};
// topic_r7_vs_r6_inverse
// lines

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
    o << content << "</body></html>";
}
