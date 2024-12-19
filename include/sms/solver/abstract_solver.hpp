#ifndef SMS_ABSTRACT_SOLVER_HPP
#define SMS_ABSTRACT_SOLVER_HPP

#include "networkit/graph/Graph.hpp"

namespace sms {

class MaxCutSolver {
public:
    explicit MaxCutSolver(const NetworKit::Graph &graph);

    virtual ~MaxCutSolver() = default;

    virtual void run() = 0;

    virtual bool applicable() = 0;

    virtual bool optimalityProven() const = 0;

    virtual double bestSolutionValue() const;

    virtual std::vector<bool> getBestSolution() const = 0;

protected:
    const NetworKit::Graph &graph_;
    double bestValue_ = 0.;

    bool hasRun_ = false;
};
} // namespace sms

#endif // SMS_ABSTRACT_SOLVER_HPP
