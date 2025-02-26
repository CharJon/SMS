#ifndef SMS_ENUMERATION_SOLVER_GRAY_HPP
#define SMS_ENUMERATION_SOLVER_GRAY_HPP

#include <cstdint>
#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/solver/abstract_solver.hpp"

namespace sms {

class EnumerationSolverGray : public MaxCutSolver {
    using partition_t = uint_fast32_t;

public:
    explicit EnumerationSolverGray(const NetworKit::Graph &g) : MaxCutSolver(g) {};

    // The graph has to be small enough to enumerate all solutions and needs to be compact.
    // Compactness is required to avoid performance issues.
    bool applicable() override;

    bool optimalityProven() const override;

    void run() override;

    std::vector<bool> getBestSolution() const override;

    bool partitionInBestSolution(NetworKit::node u) const;

private:
    partition_t bestPartition_ = 0;

    void solveNoFixed();

    double calcOffset(partition_t currentPartition, partition_t vertexToFlip) const;
};
} // namespace sms

#endif // SMS_ENUMERATION_SOLVER_GRAY_HPP
