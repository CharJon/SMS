#include "sms/graph/fast_vertex_separators.hpp"

namespace sms {

void FastVertexSeparator::run(unsigned int maxSize) {
    maxSize_ = maxSize;
    for (auto w : graph_.nodeRange()) {
        if (activeNodes_[w])
            findSeparator(w);
    }
    hasRun_ = true;
}

void FastVertexSeparator::findSeparator(NetworKit::node w) {
    std::vector<NetworKit::node> leftSide = {w};
    vertexInComponent_[w] = true;
    auto neighborhood = std::vector<NetworKit::node>();
    for (auto v : graph_.neighborRange(w)) {
        vertexInComponent_[v] = true;
        neighborhood.push_back(v);
    }

    std::vector<NetworKit::node> bestSeparator;
    std::vector<NetworKit::node> bestSeparatedSet;

    while (leftSide.size() + neighborhood.size() <= maxSize_) {
        if (neighborhood.size() <= maxSeparatorSize_) {
            bestSeparator = neighborhood;
            bestSeparatedSet = leftSide;
        }
        size_t i = findBestNeighbor(neighborhood);
        if (i == std::numeric_limits<size_t>::max())
            break;
        NetworKit::node v = neighborhood[i];
        leftSide.push_back(v);
        neighborhood[i] = neighborhood.back();
        neighborhood.pop_back();
        updateNeighborhood(neighborhood, v);
    }

    // reset vertex type
    for (auto u : leftSide) {
        vertexInComponent_[u] = false;
    }

    for (auto u : neighborhood) {
        vertexInComponent_[u] = false;
    }

    if (!bestSeparator.empty()) {
        // Update activeNodes_
        deactivateNodes(bestSeparatedSet);
        // Save separator
        separators_.push_back(bestSeparator);
        separatedSet_.push_back(bestSeparatedSet);
    }
}

size_t FastVertexSeparator::findBestNeighbor(std::vector<NetworKit::node> neighborhood) {
    size_t indexOfOptimalNeighbor = std::numeric_limits<size_t>::max();
    size_t newNeighborhoodSize = std::numeric_limits<size_t>::max();
    size_t oldNeighbourhoodSize = neighborhood.size();
    for (size_t i = 0; i < oldNeighbourhoodSize; i++) {
        NetworKit::node v = neighborhood[i];
        if (graph_.degree(v) > maxSize_) // guarantee linear runtime
            continue;
        auto calculatedNeighborhoodSize = calcNewNeighborhoodSize(static_cast<int>(oldNeighbourhoodSize), v);
        if (static_cast<size_t>(calculatedNeighborhoodSize) <= newNeighborhoodSize) {
            indexOfOptimalNeighbor = i;
            newNeighborhoodSize = calculatedNeighborhoodSize;
        }
    }
    return indexOfOptimalNeighbor;
}

int FastVertexSeparator::calcNewNeighborhoodSize(int oldNeighborhoodSize, NetworKit::node v) {
    int disjointNeighbors = 0;
    for (auto u : graph_.neighborRange(v)) {
        if (!vertexInComponent_[u])
            disjointNeighbors++;
    }
    return oldNeighborhoodSize - 1 + disjointNeighbors;
}

void FastVertexSeparator::updateNeighborhood(std::vector<NetworKit::node> &neighborhood, NetworKit::node v) {
    for (auto u : graph_.neighborRange(v)) {
        if (!vertexInComponent_[u]) {
            vertexInComponent_[u] = true;
            neighborhood.push_back(u);
        }
    }
}

void FastVertexSeparator::deactivateNodes(const std::vector<NetworKit::node> &nodesToDeactivate) {
    for (auto v : nodesToDeactivate) {
        activeNodes_[v] = false;
    }
}

} // namespace sms
