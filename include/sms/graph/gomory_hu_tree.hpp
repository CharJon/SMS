#ifndef SMS_GOMORY_HU_TREE_HPP
#define SMS_GOMORY_HU_TREE_HPP

#include <map>
#include <stdexcept>

#include "networkit/flow/EdmondsKarp.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"

namespace sms {

class GomoryHuTree {

public:
    explicit GomoryHuTree(const NetworKit::Graph &g)
        : graph_(g), gomoryHuTree_(graph_.upperNodeIdBound(), true, false) {
        bool smallestNodeFound = false;
        for (NetworKit::node i = 0; i < graph_.upperNodeIdBound(); ++i) {
            contractions_[i] = {i};
            contractedTo_[i] = i;
            if (graph_.hasNode(i) && !smallestNodeFound) {
                smallestNode_ = i;
                smallestNodeFound = true;
            }
        }
        representative_ = std::vector<NetworKit::node>(graph_.upperNodeIdBound(), smallestNode_);
        assert(graph_.hasNode(smallestNode_));
    }

    void run();

    bool hasRun() const { return run_; }

    const NetworKit::Graph &getGomoryHuTree() const;

    NetworKit::WeightedEdge minCutEdge(NetworKit::node, NetworKit::node);

    // Contract v into w
    void contractNodesInGHT(NetworKit::node, NetworKit::node);

private:
    bool run_ = false;
    const NetworKit::Graph &graph_;
    NetworKit::Graph gomoryHuTree_;
    std::vector<NetworKit::node> representative_;
    std::vector<double> flow_ = std::vector<double>(graph_.upperNodeIdBound(), 0);
    NetworKit::node smallestNode_;
    std::vector<NetworKit::node> contractedTo_ = std::vector<NetworKit::node>(graph_.upperNodeIdBound());
    std::vector<std::vector<NetworKit::node>> contractions_ =
        std::vector<std::vector<NetworKit::node>>(graph_.upperNodeIdBound());

    void constructGomoryHuTree();
};
} // namespace sms

#endif // SMS_GOMORY_HU_TREE_HPP
