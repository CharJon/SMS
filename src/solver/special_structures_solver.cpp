#include "sms/solver/special_structures_solver.hpp"

#include "sms/graph/bipartite.hpp"
#include "sms/graph/dinic.hpp"
#include "sms/graph/graphs.hpp"

namespace sms {

SpecialStructuresSolver::SpecialStructuresSolver(const NetworKit::Graph &g)
    : MaxCutSolver(g), bestPartition_(g.upperNodeIdBound(), 3) {}

bool SpecialStructuresSolver::applicable() {
    applicabilityChecked_ = true;

    if (graph_.numberOfNodes() != graph_.upperNodeIdBound())
        throw std::runtime_error("Not implemented for non compact graphs.");

    // Compute values needed later
    computeGraphProperties(graph_);

    // Complete graphs with all positive edges of equal weight are easy
    if (numPossibleEdges_ == numPositiveEdges_ && positiveEdgesEqualWeight_) {
        for (NetworKit::node i = 0; i < n_ / 2; i++) {
            bestPartition_[i] = 0;
        }
        for (NetworKit::node i = n_ / 2; i < n_; i++) {
            bestPartition_[i] = 1;
        }
        bestValue_ = positiveEdgesUnitWeight_ * static_cast<double>(numEdgesKHalfHalf_);
        applicable_ = true;
        return true;
    }

    // Graphs with only negative edges are easy to solve
    if (numPositiveEdges_ == 0) {
        auto negativeMaxCut = getNegativeEdgesOnlyMaxCut(graph_);
        if (negativeMaxCut.has_value()) {
            bestPartition_ = negativeMaxCut.value().first;
            bestValue_ = negativeMaxCut.value().second;
            applicable_ = true;
            return true;
        }
    }

    // Bipartite graphs are easy (negative edges must not cross the cut)
    if (numPositiveEdges_ <= numEdgesKHalfHalf_) {
        auto bipartition = getPositiveBipartitionNegativeSeparation(graph_);
        if (bipartition.has_value()) {
            bestPartition_ = bipartition.value();
            bestValue_ = sumPositiveEdgeWeights_;
            applicable_ = true;
            return true;
        }
    }

    if (positiveEdgesEqualWeight_ && numPositiveEdges_ >= numEdgesKHalfHalf_) {
        auto kHalfHalfPartition = getPositiveKhalfhalfPartition(graph_);

        if (kHalfHalfPartition.has_value()) {
            bestPartition_ = kHalfHalfPartition.value().first;
            bestValue_ = kHalfHalfPartition.value().second;
            applicable_ = true;
            return true;
        }
    }

    if (numPositiveEdges_ == 1) {
        auto onePositiveEdgePartition = getOnlyOnePositiveEdgeMaxCut(graph_);

        if (onePositiveEdgePartition.has_value()) {
            bestPartition_ = onePositiveEdgePartition.value().first;
            bestValue_ = onePositiveEdgePartition.value().second;
            applicable_ = true;
            return true;
        }
    }

    applicable_ = false;
    return false;
}

void SpecialStructuresSolver::run() {
    if (not applicabilityChecked_)
        applicable();

    if (not applicable_) {
        throw std::runtime_error("Solver is not applicable!");
    }

    // the actual work is done in the applicability check
    hasRun_ = true;
}

std::vector<bool> SpecialStructuresSolver::getBestSolution() const {
    if (!hasRun_)
        throw std::runtime_error("Call run() first!");
    std::vector<bool> bestSolution(graph_.upperNodeIdBound(), false);
    for (unsigned int i = 0; i < bestSolution.size(); i++) {
        bestSolution[i] = bestPartition_[i];
    }
    return bestSolution;
}

void SpecialStructuresSolver::computeGraphProperties(const NetworKit::Graph &g) {
    n_ = graph_.numberOfNodes();
    m_ = graph_.numberOfEdges();
    numPossibleEdges_ = (n_ * (n_ - 1)) / 2;
    numEdgesKHalfHalf_ = (n_ / 2) * ((n_ + 1) / 2);

    if (!g.isWeighted()) {
        positiveEdgesUnitWeight_ = 1;
        positiveEdgesEqualWeight_ = true;
        numPositiveEdges_ = g.numberOfEdges();
        numNegativeEdges_ = 0;
        sumPositiveEdgeWeights_ = static_cast<double>(g.numberOfEdges());
    } else {
        for (auto edge : g.edgeWeightRange()) {
            if (edge.weight > 0) {
                numPositiveEdges_++;
                sumPositiveEdgeWeights_ += edge.weight;
                // Check if all positive edges have same weight
                if (positiveEdgesUnitWeight_ < 0)
                    positiveEdgesUnitWeight_ = edge.weight;
                else if (edge.weight != positiveEdgesUnitWeight_)
                    positiveEdgesEqualWeight_ = false;
            } else if (edge.weight < 0) {
                numNegativeEdges_++;
            }
        }
        if (positiveEdgesEqualWeight_ && g.isWeighted()) {
            for (auto edge : g.edgeWeightRange())
                if (edge.weight > 0)
                    assert(edge.weight == positiveEdgesUnitWeight_);
        }
    }
}

std::optional<std::vector<uint8_t>>
SpecialStructuresSolver::getPositiveBipartitionNegativeSeparation(const NetworKit::Graph &g) const {
    Bipartite bipartite(g);
    bipartite.run(true);

    if (!bipartite.isBipartiteGraph())
        return {};

    return bipartite.getPartition();
}

std::optional<std::pair<std::vector<uint8_t>, double>>
SpecialStructuresSolver::getNegativeEdgesOnlyMaxCut(const NetworKit::Graph &g) const {
    assert(numPositiveEdges_ == 0);
    return std::pair<std::vector<uint8_t>, double>{std::vector<uint8_t>(g.numberOfNodes(), 0), 0.0};
}

NetworKit::Graph SpecialStructuresSolver::getPositiveEdgesOnlyComplementGraph(const NetworKit::Graph &g) const {
    if (!g.isWeighted())
        return unweightedComplementGraph(g);

    if (g.isDirected() || (g.numberOfNodes() != g.upperNodeIdBound())) {
        throw std::runtime_error("Not implemented for directed graphs or non-compact graphs.");
    }

    // This is the same code as the unweighted complement graph, except we take weight into account and ignore negative
    // weight edges
    auto complementGraph = NetworKit::Graph(g.upperNodeIdBound());

    auto lastNeighbor = std::vector<NetworKit::node>(
        g.upperNodeIdBound(),
        std::numeric_limits<NetworKit::node>::max()); // max of unsigned type here allows us to get 0 if we add 1
    for (NetworKit::node u = 0; u < g.upperNodeIdBound(); u++) {
        for (auto neighbor : g.weightNeighborRange(u)) {
            if (neighbor.second > 0) {
                for (auto v = lastNeighbor[neighbor.first] + 1; v < std::min(u, neighbor.first); v++) {
                    assert(!complementGraph.hasEdge(neighbor.first, v));
                    complementGraph.addEdge(neighbor.first, v);
                }
                lastNeighbor[neighbor.first] = u;
            }
        }
    }

    for (NetworKit::node u = 1; u < g.upperNodeIdBound(); u++) {
        for (NetworKit::node v = lastNeighbor[u] + 1; v < u; v++) {
            complementGraph.addEdge(u, v);
        }
    }

    return complementGraph;
}

std::optional<std::pair<std::vector<uint8_t>, double>>
SpecialStructuresSolver::getPositiveKhalfhalfPartition(const NetworKit::Graph &g) const {
    std::vector<uint8_t> partition(g.upperNodeIdBound(), 3);
    double partitionValue = 0;

    if ((m_ >= (numPossibleEdges_ - n_ / 2)) || (tryFullKhalfhalfSearch_ && (m_ >= numEdgesKHalfHalf_))) {
        // Graph may contain a K_{n/2, n/2}, so we get the complement graph
        NetworKit::Graph complement = getPositiveEdgesOnlyComplementGraph(g); //  getComplementGraph(graph_);

        NetworKit::ConnectedComponents components(complement);
        components.run();
        auto connComps = components.getComponents();

        std::vector<unsigned int> componentSizes(components.numberOfComponents());

        NetworKit::count compNum = 0;

        for (const std::vector<NetworKit::node> &comp : connComps) {
            if (comp.size() > (n_ + 1) / 2) {
                return {};
            }
            componentSizes[compNum] = comp.size();
            compNum++;
        }

        std::vector<std::vector<bool>> doable(components.numberOfComponents());
        for (auto &i : doable)
            i.resize((n_ + 1) / 2 + 1, false);

        // We see if the CCs in the complement graph can be divided into equal size groups
        doable[0][0] = true;
        doable[0][componentSizes[0]] = true;

        for (unsigned int i = 1; i < doable.size(); i++) {
            for (unsigned int j = 0; j < (n_ + 1) / 2 + 1; j++) {
                if (static_cast<int>(j) - static_cast<int>(componentSizes[i]) >= 0)
                    doable[i][j] = doable[i - 1][j] || doable[i - 1][j - componentSizes[i]];
                else
                    doable[i][j] = doable[i - 1][j];
            }
        }

        // Now we can check if the graph contains a K_{n/2, n/2}
        if (doable[componentSizes.size() - 1][(n_ + 1) / 2]) {
            partitionValue = positiveEdgesUnitWeight_ * static_cast<double>(numEdgesKHalfHalf_);

            NetworKit::count pos = (n_ + 1) / 2;
            for (NetworKit::count j = componentSizes.size() - 1; j > 0; j--) {
                assert(doable[j][pos]);
                if (doable[j][pos] == doable[j - 1][pos]) {
                    for (NetworKit::node v : connComps[j])
                        partition[v] = 1;
                } else {
                    for (NetworKit::node v : connComps[j])
                        partition[v] = 0;
                    pos -= connComps[j].size();
                }
            }

            if (pos == 0)
                for (NetworKit::node v : connComps[0])
                    partition[v] = 1;
            else
                for (NetworKit::node v : connComps[0])
                    partition[v] = 0;

            return std::pair<std::vector<uint8_t>, double>{partition, partitionValue};
        }
    }

    return {};
}

std::optional<std::pair<std::vector<uint8_t>, double>>
SpecialStructuresSolver::getOnlyOnePositiveEdgeMaxCut(const NetworKit::Graph &g) const {
    std::vector<uint8_t> partition(g.upperNodeIdBound(), 0);
    double partitionValue = 0;

    NetworKit::node s = 0, t = 0;
    double positiveEdgeWeight = -1;

    NetworKit::Graph gNegFlipped(g.numberOfNodes(), true);

    for (auto edge : g.edgeWeightRange()) {
        if (edge.weight > 0) {
            assert(s == 0 && t == 0);
            s = edge.v;
            t = edge.u;
            positiveEdgeWeight = edge.weight;
        } else if (edge.weight < 0) {
            gNegFlipped.addEdge(edge.u, edge.v, -1.0 * edge.weight);
        }
    }

    assert(s != 0 || t != 0);
    assert(positiveEdgeWeight > 0);

    gNegFlipped.indexEdges();
    Dinic flow(gNegFlipped);
    flow.run(s, t);

    // If the minimum cut on gNegFlipped is less than the weight of the positive edge then taking the positive edge
    // and further vertices, so we get exactly the corresponding negative edges in the cut is worth it
    if (flow.getFlowValue() < positiveEdgeWeight) {
        for (NetworKit::node v : flow.getSourceSet())
            partition[v] = 1;
        partitionValue = positiveEdgeWeight - flow.getFlowValue();
    }

    return std::pair<std::vector<uint8_t>, double>{partition, partitionValue};
}

bool SpecialStructuresSolver::optimalityProven() const {
    assert(hasRun_ && applicabilityChecked_);
    return applicable_;
}

void SpecialStructuresSolver::setTryFullKhalfhalfSearch(bool tryFullKhalfhalfSearch) {
    tryFullKhalfhalfSearch_ = tryFullKhalfhalfSearch;
}

} // namespace sms
