#include "sms/solver/heuristic_solver.hpp"

sms::HeuristicSolver::HeuristicSolver(const NetworKit::Graph &g)
    : MaxCutSolver(g),
      mqlibInstance_(mqlibEdgelist(g), static_cast<int>(g.upperNodeIdBound())),
      bestSolution_(graph_.upperNodeIdBound(), false) {
    if (g.numberOfNodes() != g.upperNodeIdBound()) {
        throw std::runtime_error("Graph is not compact!");
    }
    hasRun_ = false;
    bestValue_ = 0;
}

void sms::HeuristicSolver::run() {
    hasRun_ = true;
    mqlib::Burer2002 mqHeur(mqlibInstance_, timelimit_.count(), false, nullptr,
                            solutionValueLimit_); // constructor automatically starts heuristic

    // update best solution and best value
    const auto mqlibAssignment = mqHeur.get_best_solution().get_assignments();
    assert(mqlibAssignment.size() == graph_.numberOfNodes());
    for (size_t i = 0; i < mqlibAssignment.size(); i++) {
        assert(mqlibAssignment[i] == -1 || mqlibAssignment[i] == 1);
        bool mqLibSide = (mqlibAssignment[i] + 1) / 2;
        bestSolution_[i] = mqLibSide;
    }
    bestValue_ = mqHeur.get_best();

    solutionValues_ = mqHeur.get_past_solution_values();
    solutionTimeStamps_ = mqHeur.get_past_solution_times();
}

std::vector<bool> sms::HeuristicSolver::getBestSolution() const {
    return bestSolution_;
}

bool sms::HeuristicSolver::applicable() {
    return true;
}

bool sms::HeuristicSolver::optimalityProven() const {
    return false;
}

void sms::HeuristicSolver::setTimelimit(std::chrono::duration<double> t) {
    timelimit_ = t;
}

void sms::HeuristicSolver::setSolutionValueLimit(double value) {
    solutionValueLimit_ = value;
}
