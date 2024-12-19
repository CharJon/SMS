#ifndef SMS_FRUSTRATIONINDEXGRAPH_H
#define SMS_FRUSTRATIONINDEXGRAPH_H

#include <optional>
#include <tuple>
#include <vector>

#include "nlohmann/json.hpp"

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"

namespace sms {
/*
 * Represents a Frustration Index instance. See: https://onlinelibrary.wiley.com/doi/full/10.1002/net.21907
 * The frustration index of a signed graph (only 1 / + and -1 / - edges) is the minimal number of edges
 * that need to be removed from a graph such that the resulting graph admits a balanced bipartition.
 * A balanced bipartition of a signed graph is a partition of the nodes into two sets such that
 * all 1-edges connect nodes from the same partition and all -1-edges connect nodes from different partitions.
 */
class FrustrationIndex {
public:
    explicit FrustrationIndex(const NetworKit::Graph &g);

    int getNumberOfNodes() const { return n_; }

    int getNumberOfEdges() const { return m_; }

    double getSolutionValue(const std::vector<uint8_t> &solVector) const;

    const NetworKit::Graph &getGraph() const { return graph_; }

    NetworKit::node getOriginalNode(NetworKit::node newNode) const { return newToOriginalNode_.at(newNode); }

    NetworKit::node getNewNode(NetworKit::node originalNode) const { return originalToNewNode_.at(originalNode); }

    bool isConsistent();

    nlohmann::ordered_json getInstanceInformation();

    void printInstanceInformation(std::ostream &out);

protected:
    std::string name_;
    int n_;
    int m_;
    std::unordered_map<NetworKit::node, NetworKit::node> originalToNewNode_;
    std::vector<NetworKit::node> newToOriginalNode_;
    NetworKit::Graph graph_;
};
} // namespace sms

#endif // SMS_FRUSTRATIONINDEXGRAPH_H
