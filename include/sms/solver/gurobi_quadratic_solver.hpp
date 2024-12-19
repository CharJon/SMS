#ifndef SMS_GUROBI_QUADRATIC_SOLVER_HPP
#define SMS_GUROBI_QUADRATIC_SOLVER_HPP

#include "gurobi_c++.h"
#include <string>

#include "sms/instance/maxcut.hpp"
#include "sms/instance/qubo.hpp"
#include "sms/instance/qubo_solution.hpp"
#include "sms/solver/abstract_solver.hpp"

namespace sms {

class GurobiQuadraticSolver : public MaxCutSolver {
public:
    explicit GurobiQuadraticSolver(const NetworKit::Graph &g);

    ~GurobiQuadraticSolver() { delete vars_; };

    void run() override;

    double getRuntimeSeconds() const;

    bool optimalityProven() const override;

    double bestSolutionValue() const override;

    nlohmann::json getStats() const;

    std::vector<bool> getBestSolution() const override;

    void setTimelimit(double timelimit) { timelimit_ = timelimit; }

    void setSeed(int seed) { seed_ = seed; }

    void setStats(const std::string &stats) { stats_ = stats; }

    void setLogGurobi(const std::string &logGurobi) { logGurobi_ = logGurobi; }

    void setThreads(int threads) { threads_ = threads; }

    bool applicable() override { return true; }

private:
    int seed_ = 0;
    double timelimit_;
    int threads_ = 1;
    std::string stats_;
    std::string logGurobi_;

    GRBEnv env_;
    GRBModel model_;
    GRBVar *vars_;
};

} // namespace sms

#endif // SMS_GUROBI_QUADRATIC_SOLVER_HPP
