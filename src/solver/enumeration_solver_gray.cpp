#include "sms/solver/enumeration_solver_gray.hpp"

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/gray_code.hpp"

namespace sms {

bool EnumerationSolverGray::applicable() {
    return (graph_.numberOfNodes() < sizeof(partition_t) * CHAR_BIT)
           & (graph_.upperNodeIdBound() == graph_.numberOfNodes());
}

void EnumerationSolverGray::run() {
    if (not applicable())
        throw std::runtime_error("EnumerationSolverGray is not applicable to this graph");

    hasRun_ = true;
    if (graph_.numberOfNodes() == 0) {
        bestValue_ = 0;
        bestPartition_ = 0;
    } else {
        solveNoFixed();
    }
}

double EnumerationSolverGray::calcOffset(partition_t currentPartition, partition_t vertexToFlip) const {
    double offset = 0;
    const auto vertexPart = (currentPartition >> vertexToFlip) & 1UL;
    for (auto [neighbor, weight] : graph_.weightNeighborRange(vertexToFlip)) {
        const auto otherPart = (currentPartition >> neighbor) & 1UL;
        const auto factor = static_cast<double>(1 - (2 * (vertexPart ^ otherPart)));
        offset += factor * weight;
    }
    return offset;
}

void EnumerationSolverGray::run(const std::vector<NetworKit::node> &, const std::vector<u_int8_t> &) {
    throw std::runtime_error("Not implemented");
}

void EnumerationSolverGray::solveNoFixed() {
    partition_t currentPartition = 0;
    bestPartition_ = currentPartition;
    double currentValue = 0;
    bestValue_ = currentValue;

    for (partition_t i = 1; i < (1UL << (graph_.numberOfNodes() - 1)); i++) {
        // partition to gray code
        partition_t nextPartition = i ^ (i >> 1);
        partition_t vertexToFlip = diffPosition(currentPartition, nextPartition);
        double offset = calcOffset(currentPartition, vertexToFlip);
        currentValue = currentValue + offset;
        currentPartition = nextPartition;
        if (currentValue > bestValue_) {
            bestPartition_ = currentPartition;
            bestValue_ = currentValue;
        }
    }
}

void EnumerationSolverGray::solveFixed() {
    throw std::runtime_error("Not implemented");
}

std::vector<bool> EnumerationSolverGray::getBestSolution() const {
    assert(hasRun_);
    auto sol = std::vector<bool>(graph_.numberOfNodes(), false);
    for (auto u : graph_.nodeRange()) {
        sol[u] = (bestPartition_ >> u) & 1u;
    }
    return sol;
}

bool EnumerationSolverGray::partitionInBestSolution(NetworKit::node u) const {
    assert(hasRun_);
    return (bestPartition_ >> u) & 1u;
}

bool EnumerationSolverGray::optimalityProven() const {
    assert(hasRun_);
    return true;
}

} // namespace sms