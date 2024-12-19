#ifndef SMS_OCW_SEPERATOR_HPP
#define SMS_OCW_SEPERATOR_HPP

#include <optional>
#include <ranges>

#include "networkit/graph/Graph.hpp"
#include "tlx/container/d_ary_addressable_int_heap.hpp"

#include "sms/graph/odd_closed_walk.hpp"

namespace sms {

class OcwSeparator {
public:
    explicit OcwSeparator(const NetworKit::Graph &originalGraph);

    void updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w);

    void updateEdgeWeight(NetworKit::index idx, NetworKit::edgeweight ew);

    /*
     * Update the blacklist. Blacklisted vertices will be ignored for shortest path search.
     */
    void updateBlacklist(NetworKit::node u, bool blacklisted) {
        assert(u < originalGraph_.upperNodeIdBound());
        blacklist_[u] = blacklisted;
    }

    std::optional<OddClosedWalk> getMostViolatedOcw(NetworKit::node source);

    // path as vector of edges with bool value true if cross and false if stay
    std::vector<ParityEdge> getShortestPath(const unsigned long bestIntermediate);

    // extract simple odd cycles from a given odd closed walk
    std::vector<OddClosedWalk> simpleOddCycles(const OddClosedWalk &ocw) const;

    // extract shorter cycles from a given odd closed walk
    std::vector<OddClosedWalk> extractCycles(const OddClosedWalk &ocw) const;

    std::vector<OddClosedWalk> getViolatedOcws(NetworKit::node source);

private:
    // const
    const NetworKit::edgeweight kInfDist_ = 1.0 + 10e-4;

    // Graphs
    uint64_t const n_;
    const NetworKit::Graph &originalGraph_;
    std::vector<NetworKit::edgeweight> weights_;

    // Output
    std::vector<uint64_t> cycles_;

    // Config
    double epsilon_ = 0.0;

    // SSSP
    std::vector<NetworKit::edgeweight> distances_;
    std::vector<int> hopDistances_;
    std::vector<NetworKit::node> previous_;
    NetworKit::node currentSource_;
    std::vector<bool> blacklist_;
    std::vector<NetworKit::node> intermediateVertices_; // vertices on short path (length < 1)

    inline NetworKit::node twinId(NetworKit::node u) const { return u ^ 1; }

    // from original vertex id to left or right vertex
    inline NetworKit::node toSide(NetworKit::node u, bool side) const { return 2 * u + side; }

    inline NetworKit::node toLeft(NetworKit::node u) const { return 2 * u; }

    inline NetworKit::node toRight(NetworKit::node u) const { return 2 * u + 1; }

    inline NetworKit::node originalId(NetworKit::node u) const { return u >> 1; }

    inline bool sameSide(NetworKit::node u, NetworKit::node v) const { return !(side(u) ^ side(v)); }

    inline bool side(NetworKit::node u) const { return u & 1u; }

    void blacklistIntReachable(NetworKit::node source);

    void runDijsktra();
    OddClosedWalk intermediateVertexToOcw(NetworKit::node intermediate);
};

} // namespace sms

#endif // SMS_OCW_SEPERATOR_HPP
