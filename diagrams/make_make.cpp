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
//topic_r7_vs_r6_inverse

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
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlog.gpl \\\n";
            }
        }
        o_v << "gpls_" + topic.name + " = " + gpls + "\n\n";
        o_v << "gpls += $(gpls_" + topic.name + ")\n\n";
        o_r << "$(gpls_" + topic.name + ") &: diagram_" + topic.name + "\n\t./diagram_" + topic.name + "\n\n";
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
                gpls += '\t' + graph.ztag + '_' + topic.name + '_' + range.name + "_zlog.gpl \\\n";
            }
        }
        o_v << "gpls_" + topic.name + " = " + gpls + "\n\n";
        o_v << "gpls += $(gpls_" + topic.name + ")\n\n";
        o_r << "$(gpls_" + topic.name + ") &: diagram_" + topic.name + "\n\t./diagram_" + topic.name + "\n\n";
    }
}
