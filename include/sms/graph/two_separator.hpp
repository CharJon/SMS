#ifndef SMS_TWO_SEPERATOR_HPP
#define SMS_TWO_SEPERATOR_HPP

#include <vector>

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/components/DynConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace sms {

bool isBiconnected(NetworKit::Graph const &g);

class TwoSeparator {
public:
    NetworKit::Graph const &graph;
    std::vector<std::tuple<NetworKit::node, NetworKit::node>> separators;

    explicit TwoSeparator(const NetworKit::Graph &graph) : graph(graph) {
        assert(isBiconnected(graph));
        separators = {};
    }

    void run();

    bool hasRun() const;

private:
    bool run_ = false;
};

} // namespace sms

#endif // SMS_TWO_SEPERATOR_HPP
