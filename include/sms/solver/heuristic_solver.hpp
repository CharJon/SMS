#ifndef SMS_HEURISTIC_SOLVER_HPP
#define SMS_HEURISTIC_SOLVER_HPP

#include <chrono>

#include "mqlib/heuristics/maxcut/burer2002.h"
#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/mqlib.hpp"
#include "sms/solver/abstract_solver.hpp"

namespace sms {

class HeuristicSolver : public MaxCutSolver {
public:
    explicit HeuristicSolver(const NetworKit::Graph &g);

    // Any weighted undirected graph is a valid input
    bool applicable() override;

    void run() override;

    std::vector<bool> getBestSolution() const override;

    // Always false: Solver is heuristic
    bool optimalityProven() const override;
    ;

    void setTimelimit(std::chrono::duration<double> t);

    void setSolutionValueLimit(double value);

    const std::vector<double> &getSolutionValues() const { return solutionValues_; }

    const std::vector<double> &getSolutionTimeStamps() const { return solutionTimeStamps_; }

private:
    mqlib::MaxCutInstance mqlibInstance_;
    std::vector<bool> bestSolution_;
    std::chrono::duration<double> timelimit_{1.0};
    double solutionValueLimit_{std::numeric_limits<double>::max()};
    std::vector<double> solutionValues_;
    std::vector<double> solutionTimeStamps_;
};

} // namespace sms

#endif // SMS_HEURISTIC_SOLVER_HPP
