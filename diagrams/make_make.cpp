#include "diagram_all.hpp"

    int 
main ()
{
        auto
    o_v = std::ofstream { "make_vars.mk" };
        auto
    o_r = std::ofstream { "make_rules.mk" };
    for (auto&& topic: topics)
    {
            auto
        gpls = std::string {};
        for (auto&& range: topic.ranges)
        {
            for (auto&& graph: topic.graphs)
            {
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlin.gpl \\\n";
                if (graph.skip_zlog) continue;
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlog.gpl \\\n";
            }
        }
        o_v << "gpls_" + topic.name + " = " + gpls + "\n\n";
        o_v << "gpls += $(gpls_" + topic.name + ")\n\n";
        o_r << "$(gpls_" + topic.name + ") &: diagram_" + topic.name + "\n\t./diagram_" + topic.name + "\n\n";
        o_r << "diagram_" + topic.name + ".o: diagram_" + topic.name + ".hpp\n\n";
    }

    {
            auto const&
        topic = topic_r7_vs_r6_inverse;
            auto
        gpls = std::string {};
        for (auto&& range: topic.ranges)
        {
            for (auto&& graph: topic.graphs)
            {
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlin.gpl \\\n";
                if (graph.skip_zlog) continue;
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlog.gpl \\\n";
            }
        }
        o_v << "gpls_" + topic.name + " = " + gpls + "\n\n";
        o_v << "gpls += $(gpls_" + topic.name + ")\n\n";
        o_r << "$(gpls_" + topic.name + ") &: diagram_" + topic.name + "\n\t./diagram_" + topic.name + "\n\n";
        o_r << "diagram_" + topic.name + ".o: diagram_" + topic.name + ".hpp\n\n";
    }
}
