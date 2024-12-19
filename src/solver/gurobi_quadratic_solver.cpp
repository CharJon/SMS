#include "sms/solver/gurobi_quadratic_solver.hpp"

#include "gurobi_c++.h"
#include <fstream>

#include "nlohmann/json.hpp"

namespace sms {

GurobiQuadraticSolver::GurobiQuadraticSolver(const NetworKit::Graph &g)
    : MaxCutSolver(g), stats_(""), logGurobi_(""), env_(), model_(env_) {
    assert(graph_.checkConsistency());
}

double GurobiQuadraticSolver::getRuntimeSeconds() const {
    if (!hasRun_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    return model_.get(GRB_DoubleAttr_Runtime);
}

bool GurobiQuadraticSolver::optimalityProven() const {
    if (!hasRun_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    return model_.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
}

double GurobiQuadraticSolver::bestSolutionValue() const {
    if (!hasRun_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }

    return model_.get(GRB_DoubleAttr_ObjVal);
}

nlohmann::json GurobiQuadraticSolver::getStats() const {
    if (!hasRun_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    nlohmann::json j;

    j["cpu_solving_time"] = model_.get(GRB_DoubleAttr_Runtime);
    j["best_solution_value"] = model_.get(GRB_DoubleAttr_ObjVal);
    j["optimal_solution_found"] = model_.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
    j["timelimit_reached"] = model_.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT;
    j["bab_nodes"] = model_.get(GRB_DoubleAttr_NodeCount);
    j["obj_bound_value"] = model_.get(GRB_DoubleAttr_ObjBound);
    j["seed"] = seed_;
    j["timelimit"] = timelimit_;
    j["threads"] = threads_;

    return j;
}

void GurobiQuadraticSolver::run() {
    int numberOfNodes = static_cast<int>(graph_.numberOfNodes());

    try {
        vars_ = model_.addVars(numberOfNodes, GRB_BINARY);

        if (!logGurobi_.empty()) {
            model_.set(GRB_StringParam_LogFile, logGurobi_);
            model_.set(GRB_IntParam_LogToConsole, 0);
        }

        model_.set(GRB_IntParam_Seed, seed_);
        model_.set(GRB_DoubleParam_TimeLimit, timelimit_);
        model_.set(GRB_IntParam_Threads, threads_);
        model_.set(GRB_DoubleParam_MIPGap, 1e-6);
        model_.set(GRB_IntParam_Presolve, 2);
        // model_.set(GRB_IntParam_Disconnected, 2);

        GRBQuadExpr quadObjective(0.0);

        for (int i = 0; i < numberOfNodes; i++) {
            quadObjective.addTerm(graph_.weightedDegree(i), vars_[i], vars_[i]);
            // vars[i].set(GRB_IntAttr_BranchPriority, static_cast<int>(graph_.degree(i)));
        }

        for (int i = 0; i < numberOfNodes; i++) {
            for (int j = i + 1; j < numberOfNodes; j++) {
                auto matrixEntry = 2 * graph_.weight(i, j);
                if (matrixEntry != 0.0) {
                    quadObjective.addTerm(-matrixEntry, vars_[i], vars_[j]);
                }
            }
        }

        model_.setObjective(quadObjective, GRB_MAXIMIZE);

        model_.optimize();

        if (!stats_.empty()) {
            std::ofstream o(stats_);
            o << std::setw(4) << getStats() << std::endl;
        }

        hasRun_ = true;

    } catch (GRBException &e) {
        std::cout << "Error code = " << e.getErrorCode() << '\n' << e.getMessage() << std::endl;
        throw e;
    }
}

std::vector<bool> GurobiQuadraticSolver::getBestSolution() const {
    if (!hasRun_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }

    std::vector<bool> solution(graph_.upperNodeIdBound(), false);
    for (unsigned int i = 0; i < graph_.upperNodeIdBound(); i++)
        solution[i] = (vars_[i].get(GRB_DoubleAttr_X) > 0.5);

    return solution;
}

} // namespace sms