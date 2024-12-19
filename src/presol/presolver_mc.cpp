#include "sms/presol/presolver_mc.hpp"

#include "networkit/auxiliary/Timer.hpp"
#include "networkit/components/ConnectedComponents.hpp"

namespace sms {

void PresolverMC::run() {
    Aux::StartedTimer t;
    decompositionPhase();
    dataReductionPhase();
    t.stop();
    elapsed_ = std::chrono::duration_cast<std::chrono::milliseconds>(t.stopTime() - t.startTime());
    hasRun_ = true;
}

unsigned int PresolverMC::biconnectedComponentsDecomposition(const NetworKit::Graph &graph,
                                                             const std::vector<NetworKit::node> &toOrig) {
    assert(graph.numberOfNodes() > 0);
    assert(graph.upperNodeIdBound() == toOrig.size());
    if (graph.numberOfNodes() == 1) {
        auto u = *graph.nodeRange().begin();
        addComponent({u}, u);
        return 1;
    }

    auto bicon = BiconnectedPartition(graph);
    bicon.run();

    std::vector<NetworKit::node> dfsStack;
    dfsStack.push_back(0); // 0 is guaranteed to be a component

    std::vector<bool> visited(bicon.tree().upperNodeIdBound(), false);
    visited[0] = true;

    auto firstComponent = bicon.components()[0];
    // transform to original node ids
    for (auto &v : firstComponent)
        v = toOrig[v];
    addComponent(firstComponent, firstComponent[0]); // any root is fine

    while (!dfsStack.empty()) {
        auto curNode = dfsStack.back();
        dfsStack.pop_back();
        for (auto neighbor : bicon.tree().neighborRange(curNode)) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                dfsStack.push_back(neighbor);

                if (bicon.isArticulationPoint(curNode)) { // articulation points have component neighbors
                    assert(!bicon.isArticulationPoint(neighbor));
                    visited[neighbor] = true;
                    dfsStack.push_back(neighbor);
                    std::vector<NetworKit::node> component = bicon.components()[neighbor];
                    for (auto &v : component)
                        v = toOrig[v];
                    addComponent(component, toOrig[bicon.correspondingArticulationVertex(curNode)]);
                }
            }
        }
    }

    return static_cast<int>(bicon.numComponents());
}

void PresolverMC::decompositionPhase() {
    auto connected = NetworKit::ConnectedComponents(graph_);
    connected.run();

    auto allComponents = connected.getComponents();
    for (unsigned int i = 0; i < connected.numberOfComponents(); i++) {
        // decompose each connected component individually
        auto inducedSubgraph = inducedSubgraphCompact(graph_, allComponents[i]);
        biconnectedComponentsDecomposition(inducedSubgraph.compactGraph, inducedSubgraph.compact2orig);
    }
    assert(subgraphs_.size() == subgraphVertexToOriginalVertex_.size());
    assert(roots_.size() == subgraphs_.size());
}

void PresolverMC::dataReductionPhase() {
    for (int i = 0; i < std::ssize(reducers_); i++) {
        auto &curReducer = reducers_[i];
        curReducer.run();
    }
}

std::vector<bool> PresolverMC::recoverSolution() {
    // Recover solutions for subgraphs via data reduction history
    // Merge solutions through the subgraphs roots
    std::vector<Partition> fullSolution(graph_.upperNodeIdBound(), Partition::kUNASSIGNED);

    for (unsigned int i = 0; i < subgraphs_.size(); i++) {
        // Reverse reduction history here to get solution for full redSubgraph
        auto &reducedSolution = solutions_[i];
        auto &reducer = reducers_[i];
        auto subgraphSolution = reducer.recoverSolution(reducedSolution);

        // Perform solution merge, based on roots
        auto &sub2orig = subgraphVertexToOriginalVertex_[i];
        auto root = roots_[i];
        if (fullSolution[sub2orig[root]] == Partition::kUNASSIGNED) {
            fullSolution[sub2orig[root]] = Partition::kZERO;
        }
        bool flip = fullSolution[sub2orig[root]] != subgraphSolution[root];
        for (auto u : subgraphs_[i].nodeRange()) {
            auto origNode = sub2orig[u];
            assert(subgraphSolution[u] == Partition::kZERO || subgraphSolution[u] == Partition::kONE);
            fullSolution[origNode] = static_cast<sms::Partition>(flip ^ subgraphSolution[u]);
        }
    }
    std::vector<bool> fullSolutionBool(graph_.upperNodeIdBound(), false);
    for (unsigned int i = 0; i < fullSolution.size(); i++) {
        fullSolutionBool[i] = fullSolution[i] == Partition::kONE;
    }
    return fullSolutionBool;
}

void PresolverMC::addComponent(const std::vector<NetworKit::node> &component, NetworKit::node root) {
    assert(component.size() > 0);
    assert(std::find(component.begin(), component.end(), root) != component.end());
    assert(std::all_of(component.begin(), component.end(), [this](NetworKit::node u) { return graph_.hasNode(u); }));

    auto [inducedSubgraph, sub2orig, orig2sub] = sms::inducedSubgraphCompact(graph_, component);
    reducers_.emplace_back(inducedSubgraph, emphasis_); // before move!
    if (triangleEmphasis_ < 3) {
        reducers_.back().turnOffNewTriangles();
        if (triangleEmphasis_ < 1) {
            reducers_.back().turnOffTriangles();
        }
    }
    if (oldPresolve_)
        reducers_.back().turnOffNewPresolve();

    subgraphs_.push_back(std::move(inducedSubgraph));
    subgraphVertexToOriginalVertex_.push_back(std::move(sub2orig));
    roots_.push_back(orig2sub[root]);
}

nlohmann::ordered_json PresolverMC::getStats() const {
    nlohmann::ordered_json j;
    for (unsigned int i = 0; i < reducers_.size(); i++) {
        j["reducer " + std::to_string(i)] = reducers_[i].getStats();
    }
    return j;
}

nlohmann::ordered_json PresolverMC::getCompactStats() const {
    nlohmann::ordered_json j;
    if (reducers_.empty())
        return j;
    j = reducers_[0].getStats();
    for (unsigned int i = 1; i < reducers_.size(); i++) {
        // merge all entries of json
        auto otherStats = reducers_[i].getStats();
        for (auto it = otherStats.begin(); it != otherStats.end(); ++it) {
            if (it.key() != "emphasis")
                j[it.key()] = j[it.key()].get<unsigned int>() + it.value().get<unsigned int>();
        }
    }
    return j;
}

void PresolverMC::setOldPresolve() {
    oldPresolve_ = true;
}

} // namespace sms