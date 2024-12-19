#ifndef SMS_ENUMERATION_SOLVER_HPP
#define SMS_ENUMERATION_SOLVER_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/solver/abstract_solver.hpp"

namespace sms {

/*
 * A class to solve max cut via full enumeration of all solutions.
 * Only works for instances with less than 32 vertices.
 */
class EnumerationSolver : public MaxCutSolver {
    using partition_t = uint32_t;

public:
    explicit EnumerationSolver(const NetworKit::Graph &g) : MaxCutSolver(g) {};

    bool applicable() override;

    bool optimalityProven() const override;

    void run() override;

    /**
     * @param fixedNodes Vector of (0-indexed) node ids which should be fixed to a certain side
     * @param fixedValues Vector of fixed assignment, one entry for each node in @fixedNodes
     */
    void run(const std::vector<NetworKit::node> &fixedNodes, const std::vector<u_int8_t> &fixedValues);

    /**
     * @param fixedNodes Vector of (0-indexed) node ids which should be fixed to a certain side
     * @param fixedValues Vector of fixed assignment, one entry for each node in @fixedNodes
     */
    void fix(const std::vector<NetworKit::node> &fixedNodes, const std::vector<u_int8_t> &fixedValues);

    std::vector<bool> getBestSolution() const override;

    bool partitionOf(NetworKit::node u) const { return (bestPartition_ >> u) & 1UL; }

    NetworKit::edgeweight getCutValue(partition_t) const;

    NetworKit::edgeweight getCutValueLambda(partition_t cut) const;

private:
    partition_t bestPartition_ = 0;

    partition_t fixedValues_ = 0;
    partition_t fixedMask_ = 0;

private:
    void solveNoFixed();

    void solveFixed();
};

} // namespace sms

#endif // SMS_ENUMERATION_SOLVER_HPP
