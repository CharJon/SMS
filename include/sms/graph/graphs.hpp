#ifndef SMS_GRAPHS_HPP
#define SMS_GRAPHS_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/math.hpp"

namespace sms {

struct CompactGraph {
    NetworKit::Graph compactGraph;
    std::vector<NetworKit::node> compact2orig;
    std::unordered_map<NetworKit::node, NetworKit::node> orig2compact;
};

/*
 * Calculates a divisor for all edge weights.
 * For integers, it's the greatest common divisor of all edge weights.
 * For double only a simple heuristic is used.
 */
NetworKit::edgeweight edgeWeightDivisor(const NetworKit::Graph &graph);

/*
 * If all weights are integer:
 * Degrees are even and all weights are integer and odd, return 2.0 else 1.0
 */
NetworKit::edgeweight degreeBasedScaling(const NetworKit::Graph &graph);

/*
 * Returns the node with the maximum degree in the graph and its degree.
 * For graphs without nodes this has undefined behavior.
 */
std::tuple<NetworKit::node, NetworKit::count> maxDegreeNode(const NetworKit::Graph &g);

/*
 * Get a compact induced subgraph.
 * Only returns subgraph. New node IDs are continuous and start at 0 for the first node in @nodes etc.
 */
CompactGraph inducedSubgraphCompact(const NetworKit::Graph &g, const std::vector<NetworKit::node> &nodes);

/*
 * Get graph with compact vertex IDs.
 * New node IDs are continuous and start at 0 for the first node in @nodes etc.
 * If seed > 0, the nodes are shuffled with the given seed.
 */
CompactGraph compactGraph(const NetworKit::Graph &graph, int seed = 0);

/*
 * Get the complement graph of an input graph (complement is unweighted).
 */
NetworKit::Graph unweightedComplementGraph(const NetworKit::Graph &graph);

NetworKit::edgeweight neighborhoodAlpha(const NetworKit::Graph &g, NetworKit::node u, NetworKit::node v,
                                        bool ignoreEachOther, std::vector<NetworKit::edgeweight> &marks);

NetworKit::edgeweight neighborhoodAlphaMarksSet(const NetworKit::Graph &g, NetworKit::node u, NetworKit::node v,
                                                bool ignoreEachOther, std::vector<NetworKit::edgeweight> &marks);

NetworKit::edgeweight getCutValue(const NetworKit::Graph &g, const std::vector<bool> &solution);

/*
 * Returns true if all vertices have integer weighted degree, false otherwise.
 */
bool integerWeightedDegree(const NetworKit::Graph &g);

double integerScalar(const NetworKit::Graph &g, double maxScalar = maxSafeInteger<double>(), double eps = 10e-6);

NetworKit::Graph scaleToInt(NetworKit::Graph g, double scalar);

} // namespace sms

#endif // SMS_GRAPHS_HPP
