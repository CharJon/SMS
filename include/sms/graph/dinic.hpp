#ifndef SMS_DINIC_HPP
#define SMS_DINIC_HPP

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/math.hpp"

namespace sms {

// Dinic MaxFlow-algorithm for undirected graphs
class Dinic {
public:
    explicit Dinic(const NetworKit::Graph &g) : graph_(g), flow_(graph_.upperEdgeIdBound(), 0.0) {
        assert(graph_.checkConsistency());
    }

    void run(NetworKit::node source, NetworKit::node sink);

    bool hasRun() const { return hasRun_; }

    double getFlowValue() const { return maxFlowValue_; }

    std::vector<NetworKit::node> getSourceSet();

private:
    const NetworKit::Graph &graph_;

    NetworKit::node source_ = NetworKit::none;
    NetworKit::node sink_ = NetworKit::none;
    bool hasRun_ = false;

    double maxFlowValue_;

    int sourceLevel_ = -1;

    std::vector<int> nodeLevel_;

    std::vector<NetworKit::edgeweight> flow_; // positive means from smaller to greater, negative other way around

    bool pathFinding();
    bool rankNodes();
    void resetFlow() { std::fill(flow_.begin(), flow_.end(), 0.0); }

    // Returns current flow from u to v
    NetworKit::edgeweight flow(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edgeId) {
        // In an undirected network, edges only carry flow in one direction

        if (flow_[edgeId] > 0 && u < v)
            return flow_[edgeId];

        if (flow_[edgeId] < 0 && u > v)
            return std::fabs(flow_[edgeId]);

        return 0;
    }

    // Returns residual capacity of an edge
    NetworKit::edgeweight residualCapacity(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight capacity,
                                           NetworKit::edgeid edgeId) {
        // if flow flows the same direction, we can push the remaining capacity
        // else we can push capacity + flow in other direction
        if (u < v)
            return capacity - flow_[edgeId];
        else
            return capacity + flow_[edgeId];
    }
};

} // namespace sms

#endif // SMS_DINIC_HPP
