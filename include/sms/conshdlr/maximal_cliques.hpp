#ifndef SMS_MAXIMAL_CLIQUES_HPP
#define SMS_MAXIMAL_CLIQUES_HPP

#include "networkit/graph/Graph.hpp"
#include <networkit/base/Algorithm.hpp>

namespace sms {

/**
 * Algorithm for listing all maximal cliques.
 *
 * The implementation is based on the "hybrid" algorithm described in
 *
 * Eppstein, D., & Strash, D. (2011).
 * Listing All Maximal Cliques in Large Sparse Real-World Graphs.
 * In P. M. Pardalos & S. Rebennack (Eds.),
 * Experimental Algorithms (pp. 364 - 375). Springer Berlin Heidelberg.
 * Retrieved from http://link.springer.com/chapter/10.1007/978-3-642-20662-7_31
 *
 * The running time of this algorithm should be in @f$O(d^2 \cdot n \cdot 3^{d/3})@f$
 * where @f$d@f$ is the degeneracy of the graph, i.e., the maximum core number.
 * The running time in practice depends on the structure of the graph. In
 * particular for complex networks it is usually quite fast, even graphs with
 * millions of edges can usually be processed in less than a minute.
 */
class MaximalCliques final : public NetworKit::Algorithm {

public:
    /**
     * Construct the maximal cliques algorithm with the given graph and a callback.
     *
     * The callback is called once for each found clique with a reference to the clique.
     * Note that the reference is to an internal object, the callback should not assume that
     * this reference is still valid after it returned.
     *
     * @param G The graph to list cliques for
     * @param callback The callback to call for each clique.
     */
    MaximalCliques(const NetworKit::Graph &G, std::function<bool(const std::vector<NetworKit::node> &)> callback);

    /**
     * Execute the maximal clique listing algorithm.
     */
    void run() override;

private:
    const NetworKit::Graph *G;

    std::function<bool(const std::vector<NetworKit::node> &)> callback;
};
} // namespace sms

#endif // SMS_MAXIMAL_CLIQUES_HPP
