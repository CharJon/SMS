#ifndef SMS_MST_HEURISTIC_HPP
#define SMS_MST_HEURISTIC_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"
#include "sms/graph/odd_closed_walk.hpp"
#include "sms/instance/mc_solution.hpp"

namespace sms {

class MSTHeuristic {
private:
    NetworKit::Graph const *const originalGraph_;
    NetworKit::Graph lpWeightedGraph_;
    NetworKit::Graph tree_;

    std::queue<NetworKit::node> queue_;
    std::vector<NetworKit::node> pred_;

    std::vector<bool> flipSolution(std::vector<bool>, NetworKit::node, NetworKit::node);

public:
    explicit MSTHeuristic(NetworKit::Graph const *const originalGraph)
        : originalGraph_(originalGraph),
          lpWeightedGraph_(*originalGraph, true, false),
          tree_(originalGraph->upperNodeIdBound(), true, false),
          pred_(tree_.numberOfNodes(), -1) {}

    void updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w);

    void computeSpanningTree();

    // get primal solution from spanning tree, starting vertex s is in partition 0
    std::vector<bool> getPrimalSolution(NetworKit::node s);

    void subtreeFlip(std::vector<bool> &);

    std::vector<NetworKit::node> getSTPath(NetworKit::node u, NetworKit::node v);
};

} // namespace sms

#endif // SMS_MST_HEURISTIC_HPP
