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

    void run(const std::vector<NetworKit::node> & /*fixedNodes*/, const std::vector<u_int8_t> & /*fixedValues*/);

    std::vector<bool> getBestSolution() const override;

    bool partitionInBestSolution(NetworKit::node u) const;

private:
    partition_t bestPartition_ = 0;

    void solveNoFixed();

    void solveFixed();

    double calcOffset(partition_t currentPartition, partition_t vertexToFlip) const;
};
} // namespace sms

#endif // SMS_ENUMERATION_SOLVER_GRAY_HPP
