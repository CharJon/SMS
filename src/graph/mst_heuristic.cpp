#include "sms/graph/mst_heuristic.hpp"

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/union_find.hpp"
#include "sms/graph/odd_closed_walk.hpp"

namespace sms {

void MSTHeuristic::updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
    assert(originalGraph_->hasEdge(u, v));
    lpWeightedGraph_.setWeight(u, v, w);
}

void MSTHeuristic::computeSpanningTree() {
    tree_.removeAllEdges();
    std::vector<NetworKit::WeightedEdge> lpWeightedEdges;
    lpWeightedEdges.reserve(lpWeightedGraph_.numberOfEdges());

    lpWeightedGraph_.forEdges([&lpWeightedEdges](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        lpWeightedEdges.emplace_back(u, v, weight);
    });
    std::sort(lpWeightedEdges.begin(), lpWeightedEdges.end(),
              [](auto &e1, auto &e2) { return std::abs(e1.weight - 0.5) > std::abs(e2.weight - 0.5); });

    UnionFind unionFind(originalGraph_->upperNodeIdBound());

    for (auto edge : lpWeightedEdges) {
        if (unionFind.find(edge.u) != unionFind.find(edge.v)) {
            tree_.addEdge(edge.u, edge.v, edge.weight);
            unionFind.unite(edge.u, edge.v);
        }
    }

    assert(tree_.numberOfNodes() == originalGraph_->numberOfNodes());
    assert(tree_.numberOfEdges() == (originalGraph_->numberOfNodes() - 1));
}

std::vector<bool> MSTHeuristic::getPrimalSolution(NetworKit::node s) {
    std::vector<bool> solution(originalGraph_->upperNodeIdBound(), false);
    std::vector<bool> visited(tree_.numberOfNodes(), false);

    NetworKit::node current = s;
    queue_ = {};
    visited[current] = true;
    queue_.push(current);

    while (!queue_.empty()) {
        current = queue_.front();
        queue_.pop();
        for (auto [v, weight] : tree_.weightNeighborRange(current)) {
            if (!visited[v]) {
                queue_.push(v);
                visited[v] = true;
                solution[v] = solution[current] ^ (weight > 0.5);
            }
        }
    }
    return solution;
}

std::vector<NetworKit::node> MSTHeuristic::getSTPath(NetworKit::node u, NetworKit::node v) {
    const auto kUnassigned = std::numeric_limits<NetworKit::node>::max();
    std::fill(pred_.begin(), pred_.end(), kUnassigned);
    queue_ = {};

    queue_.push(u);
    pred_[u] = u;

    while (!queue_.empty() && pred_[v] == kUnassigned) {
        NetworKit::node current = queue_.front();
        queue_.pop();

        for (auto w : tree_.neighborRange(current)) {
            if (pred_[w] == kUnassigned) {
                queue_.push(w);
                pred_[w] = current;
            }
        }
    }

    std::vector<NetworKit::node> path;
    NetworKit::node current = v;
    while (current != u) {
        path.push_back(current);
        current = pred_[current];
    }

    path.push_back(u);

    return path;
}

void MSTHeuristic::subtreeFlip(std::vector<bool> &solution) {
    auto n = lpWeightedGraph_.numberOfNodes();
    std::vector<NetworKit::WeightedEdge> treeEdges;
    treeEdges.reserve(tree_.numberOfEdges());

    tree_.forEdges([&treeEdges](const NetworKit::node u, const NetworKit::node v, const NetworKit::edgeweight weight) {
        treeEdges.emplace_back(u, v, weight);
    });

    std::partial_sort(treeEdges.begin(), treeEdges.begin() + n / 2, treeEdges.end(),
                      [](auto e1, auto e2) { return std::abs(e1.weight - 0.5) > std::abs(e2.weight - 0.5); });

    auto currentSolutionValue = getCutValue(*originalGraph_, solution);
    std::vector<bool> currentSolution;
    currentSolution.reserve(originalGraph_->upperNodeIdBound());
    for (unsigned i = 0; i < n / 2; i++) {
        auto edge = treeEdges[i];
        currentSolution = flipSolution(solution, edge.u, edge.v);
        auto flipSolValue = getCutValue(*originalGraph_, currentSolution);
        if (flipSolValue > currentSolutionValue) {
            solution = currentSolution;
            currentSolutionValue = flipSolValue;
        }
    }
}

std::vector<bool> MSTHeuristic::flipSolution(std::vector<bool> solution, NetworKit::node u, NetworKit::node v) {
    std::vector<bool> visited(tree_.numberOfNodes(), false);
    queue_ = {};

    visited[u] = true;
    NetworKit::node current = v;
    queue_.push(current);
    visited[current] = true;
    solution[current] = !solution[current];

    while (!queue_.empty()) {
        current = queue_.front();
        queue_.pop();
        for (auto neighbor : tree_.neighborRange(current)) {
            if (!visited[neighbor]) {
                queue_.push(neighbor);
                visited[neighbor] = true;
                solution[neighbor] = !solution[neighbor];
            }
        }
    }
    return solution;
}

} // namespace sms