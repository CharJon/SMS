#include "sms/graph/two_separator.hpp"

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/DynConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace sms {

bool TwoSeparator::hasRun() const {
    return run_;
}

void TwoSeparator::run() {
    NetworKit::Graph gCopy(graph);

    NetworKit::DynConnectedComponents generator(gCopy);
    generator.run();

    for (auto u : graph.nodeRange()) {
        assert(gCopy.numberOfNodes() == graph.numberOfNodes());
        assert(gCopy.numberOfEdges() == graph.numberOfEdges());

        for (auto n : graph.neighborRange(u)) {
            gCopy.removeEdge(u, n);
            generator.update({NetworKit::GraphEvent::EDGE_REMOVAL, u, n});
        }

        for (auto v : graph.nodeRange()) {
            if (u < v) {
                for (auto n : graph.neighborRange(v)) {
                    if (n != u) {
                        gCopy.removeEdge(v, n);
                        generator.update({NetworKit::GraphEvent::EDGE_REMOVAL, v, n});
                    }
                }

                if (generator.numberOfComponents() > 3)
                    separators.emplace_back(u, v);

                for (auto n : graph.neighborRange(v)) {
                    if (n != u) {
                        gCopy.addEdge(v, n);
                        generator.update({NetworKit::GraphEvent::EDGE_ADDITION, v, n});
                    }
                }
            }
        }

        for (auto n : graph.neighborRange(u)) {
            gCopy.addEdge(u, n);
            generator.update({NetworKit::GraphEvent::EDGE_ADDITION, u, n});
        }
    }

    run_ = true;
}

bool isBiconnected(const NetworKit::Graph &g) {
    NetworKit::BiconnectedComponents generator(g);
    generator.run();
    return generator.getComponents().size() <= 1;
}

} // namespace sms
