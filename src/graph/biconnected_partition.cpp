#include "sms/graph//biconnected_partition.hpp"

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace sms {

void BiconnectedPartition::run() {
    NetworKit::BiconnectedComponents generator(graph_);
    generator.run();

    components_ = std::vector<std::vector<NetworKit::node>>(generator.numberOfComponents());

    // component tree has one node for each component
    tree_ = NetworKit::Graph(components_.size(), false, false);

    for (auto u : graph_.nodeRange()) {
        const auto &curComps = generator.getComponentsOfNode(u);
        assert(!curComps.empty());

        // check if separating vertex
        if (curComps.size() > 1) {
            auto newTreeNode = tree_.addNode(); // add node to tree
            articulationPoints_.push_back(u);
            assert(articulationPoints_[newTreeNode - components_.size()] == u);
            for (auto componentIdx : curComps) {
                components_[componentIdx].push_back(u);
                tree_.addEdge(componentIdx, newTreeNode);
            }
        } else {
            components_[*curComps.begin()].push_back(u);
        }
    }
}

} // namespace sms