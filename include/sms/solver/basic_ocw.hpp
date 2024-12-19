#ifndef SMS_BASIC_OCW_HPP
#define SMS_BASIC_OCW_HPP

#include <utility>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/scip.hpp"
#include "sms/auxiliary/scip_exception.hpp"
#include "sms/instance/maxcut.hpp"
#include "sms/instance/mc_solution.hpp"
#include "sms/solver/abstract_solver.hpp"

namespace sms {

class BasicOcwSolver : public MaxCutSolver {
public:
    explicit BasicOcwSolver(const NetworKit::Graph &graph);

    ~BasicOcwSolver();

    void run() override;

    void addWarmStartSolution(const std::vector<bool> &solution);

    double bestSolutionValue() const override;

    std::vector<bool> getBestSolution() const override;

    void setIntegral(bool integral) { integral_ = integral; }

    bool getIntegral() const { return integral_; }

    bool optimalityProven() const override;

    nlohmann::ordered_json getStats() const;

    void setTimelimit(double timelimit) { timelimit_ = timelimit; }

    void setSeed(int seed) { seed_ = seed; }

    void setNodelimit(int nodelimit) { nodelimit_ = nodelimit; }

    void setNodeVars(int nodeVars) { nodeVars_ = nodeVars; }

    void setParamFile(std::string paramFile) { paramFile_ = std::move(paramFile); }

    bool applicable() override { return true; }

private:
    SCIP *scip_;
    int seed_;
    double timelimit_;
    int nodelimit_ = -1;
    std::string paramFile_;
    bool nodeVars_ = true;

    std::vector<std::vector<bool>> warmStartSolutions_;

    bool integral_;

    std::vector<uint8_t> bestPartition_;

    void setScipParams();

    void addWarmStartSolutionsToScip();

    void addWarmStartSolutionsEdgesOnlyToScip();

    std::vector<bool> getBestSolutionEdgesOnly() const;

    std::vector<bool> solveRooted(const std::string &statsScip = "", const std::string &logScip = "",
                                  const std::string &history = "");

    std::vector<bool> solveEdgesOnly(const std::string &statsScip = "", const std::string &logScip = "",
                                     const std::string &history = "");
};

} // namespace sms
#endif // SMS_BASIC_OCW_HPP
