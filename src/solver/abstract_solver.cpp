#include "sms/solver/abstract_solver.hpp"

#include "sms/graph/graphs.hpp"

namespace sms {

MaxCutSolver::MaxCutSolver(const NetworKit::Graph &graph) : graph_(graph) {};

double MaxCutSolver::bestSolutionValue() const {
    assert(hasRun_);
    return bestValue_;
}
} // namespace sms