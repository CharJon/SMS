#ifndef SMS_BICONNECTED_PARTITION_HPP
#define SMS_BICONNECTED_PARTITION_HPP

#include <map>
#include <unordered_map>
#include <vector>

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace sms {

/*
 * Stores the biconnected components of a graph as a tree.
 * The tree is a graph where each node represents a biconnected component.
 */
class BiconnectedPartition {
public:
    explicit BiconnectedPartition(const NetworKit::Graph &graph) : graph_(graph) { tree_ = NetworKit::Graph(0); }

    void run();

    const std::vector<std::vector<NetworKit::node>> &components() const { return components_; }

    size_t numComponents() const { return components_.size(); }

    const std::vector<NetworKit::node> &articulationPoints() const { return articulationPoints_; }

    size_t numArticulationPoints() const { return articulationPoints_.size(); }

    const NetworKit::Graph &tree() const { return tree_; }

    /*
     * Given a tree node, does this tree node correspond to an articulation vertex in the original graph.
     * If false the tree node corresponds to a biconnected component.
     */
    bool isArticulationPoint(NetworKit::node u) {
        assert(tree_.hasNode(u));
        return u >= components_.size();
    }

    /*
     * Given a tree node, get the corresponding articulation vertex
     */
    NetworKit::node correspondingArticulationVertex(NetworKit::node u) {
        assert(isArticulationPoint(u));
        return articulationPoints_[u - components_.size()];
    }

private:
    NetworKit::Graph const &graph_;
    std::vector<std::vector<NetworKit::node>> components_;
    std::vector<NetworKit::node> articulationPoints_;
    // The tree storing the bicon tree
    // Two types of nodes:
    //   Node IDs [0, #components-1] are biconnected components
    //   Node IDs [#components, #components+#articulationPoints-1] are articulation points
    NetworKit::Graph tree_;
};

} // namespace sms

#endif // SMS_BICONNECTED_PARTITION_HPP
