#include "sms/solver/basic_ocw.hpp"

#include "sms/auxiliary/scip.hpp"
#include "sms/auxiliary/scip_exception.hpp"
#include "sms/branch/branchrule_degree.hpp"
#include "sms/branch/branchrule_nm.hpp"
#include "sms/conshdlr/conshdlr_cliques.hpp"
#include "sms/conshdlr/conshdlr_holes.hpp"
#include "sms/conshdlr/conshdlr_oscw.hpp"
#include "sms/conshdlr/conshdlr_triangles.hpp"
#include "sms/disp/disp.hpp"
#include "sms/eventhdlr/eventhdlr_history.hpp"
#include "sms/eventhdlr/eventhdlr_root.hpp"
#include "sms/graph/graphs.hpp"
#include "sms/instance/mc_solution.hpp"
#include "sms/pheur/pheur_kl.hpp"
#include "sms/pheur/pheur_mqlib.hpp"
#include "sms/pheur/pheur_mst.hpp"
#include "sms/probdata/mc_edges_only.hpp"
#include "sms/probdata/mc_rooted.hpp"
#include <sms/pheur/pheur_burer.hpp>

namespace sms {

BasicOcwSolver::BasicOcwSolver(const NetworKit::Graph &graph) : MaxCutSolver(graph), scip_(nullptr), paramFile_() {
    assert(graph_.numberOfNodes() > 0);
}

BasicOcwSolver::~BasicOcwSolver() {
    SCIP_CALL_EXC_NO_THROW(SCIPfree(&scip_))
    assert(scip_ == nullptr);
#ifndef NDEBUG
    BMScheckEmptyMemory();
#endif
}

void BasicOcwSolver::run() {
    if (nodeVars_) {
        solveRooted();
    } else {
        solveEdgesOnly();
    }
}

std::vector<bool> BasicOcwSolver::solveEdgesOnly(const std::string &statsScip, const std::string &logScip,
                                                 const std::string &history) {

    std::cout << "Started solver with no node-vars" << std::endl;

    hasRun_ = false;
    buildEdgesOnlyMcModel(&scip_, &graph_);
    SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, true);
    SCIPsetHeuristics(scip_, SCIP_PARAMSETTING_OFF, true);
    SCIPsetSeparating(scip_, SCIP_PARAMSETTING_OFF, true);
    if (integral_) {
        setObjectiveIntegral(scip_);
    }

    if (!logScip.empty()) {
        SCIPsetMessagehdlrLogfile(scip_, logScip.c_str());
    }

    EventhdlrHistory *eventHandler = nullptr;
    if (!history.empty()) {
        eventHandler = new EventhdlrHistory(scip_);
        SCIP_CALL_EXC(SCIPincludeObjEventhdlr(scip_, eventHandler, TRUE))
    }

    if (!statsScip.empty()) {
        std::filesystem::path p(statsScip);
        auto extension = p.extension();
        p.replace_extension(".root" + extension.string());
        auto eventHandlerRoot = new EventhdlrRoot(scip_, p.string());
        SCIP_CALL_EXC(SCIPincludeObjEventhdlr(scip_, eventHandlerRoot, TRUE))
    }

    // Add odd cycle ocwHandler
    auto *ocwFastHandler = new ConshdlrOscw(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, ocwFastHandler, TRUE));
    // Add a constraint for the odd cycle ocwHandler, to trigger callbacks to it
    SCIP_CONS *cons;
    SCIP_CALL_EXC(SCIPcreateConsBasicOddCyclesFast(scip_, &cons, "Odd Cycles", &graph_))
    SCIP_CALL_EXC(SCIPaddCons(scip_, cons))
    SCIP_CALL_EXC(SCIPreleaseCons(scip_, &cons))

    // Add triangle handler
    auto *triangleHandler = new ConshdlrTriangles(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, triangleHandler, TRUE))

    // Add hole handler
    auto *holeHandler = new ConshdlrHoles(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, holeHandler, TRUE))

    auto *cliqueHandler = new ConshdlrCliques(scip_, graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, cliqueHandler, TRUE))

    addWarmStartSolutionsEdgesOnlyToScip();
    setScipParams();

    // Run the solver
    SCIP_CALL_EXC(SCIPsolve(scip_));

    if (!statsScip.empty()) {
        scipWriteStatistics(scip_, statsScip);
    }

    if (!history.empty() && eventHandler != nullptr) {
        auto eventHandlerLog = std::ofstream(history);
        eventHandler->to_json(eventHandlerLog);
        eventHandlerLog.close();
    }

    hasRun_ = true;
    std::cout << "Solver finished." << std::endl;

    return getBestSolutionEdgesOnly();
}

std::vector<bool> BasicOcwSolver::solveRooted(const std::string &statsScip, const std::string &logScip,
                                              const std::string &history) {
    hasRun_ = false;

    auto [maxDegNode, maxDeg] = maxDegreeNode(graph_);

    // Build model
    NetworKit::node root = maxDegNode;
    std::cout << "Picked " << root << " as root node. Degree is " << maxDeg << std::endl;
    buildRootMcModel(&scip_, &graph_, root, false);
    SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, true);
    SCIPsetHeuristics(scip_, SCIP_PARAMSETTING_OFF, true);
    SCIPsetSeparating(scip_, SCIP_PARAMSETTING_OFF, true);
    if (integral_) {
        setObjectiveIntegral(scip_);
    }

    if (!logScip.empty()) {
        SCIPsetMessagehdlrLogfile(scip_, logScip.c_str());
    }

    EventhdlrHistory *eventHandler = nullptr;
    if (!history.empty()) {
        eventHandler = new EventhdlrHistory(scip_);
        SCIP_CALL_EXC(SCIPincludeObjEventhdlr(scip_, eventHandler, TRUE))
    }

    if (!statsScip.empty()) {
        std::filesystem::path p(statsScip);
        auto extension = p.extension();
        p.replace_extension(".root" + extension.string());
        auto eventHandlerRoot = new EventhdlrRoot(scip_, p.string());
        SCIP_CALL_EXC(SCIPincludeObjEventhdlr(scip_, eventHandlerRoot, TRUE))
    }

    // Add odd cycle ocwHandler
    auto *ocwFastHandler = new ConshdlrOscw(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, ocwFastHandler, TRUE));
    // Add a constraint for the odd cycle ocwHandler, to trigger callbacks to it
    SCIP_CONS *cons;
    SCIP_CALL_EXC(SCIPcreateConsBasicOddCyclesFast(scip_, &cons, "Odd Cycles", &graph_))
    SCIP_CALL_EXC(SCIPaddCons(scip_, cons))
    SCIP_CALL_EXC(SCIPreleaseCons(scip_, &cons))

    // Add triangle handler
    auto *triangleHandler = new ConshdlrTriangles(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, triangleHandler, TRUE))

    // Add hole handler
    auto *holeHandler = new ConshdlrHoles(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, holeHandler, TRUE))

    auto *cliqueHandler = new ConshdlrCliques(scip_, graph_);
    SCIP_CALL_EXC(SCIPincludeObjConshdlr(scip_, cliqueHandler, TRUE))

    // Add primal heuristics
    // auto pHeurMqLib = new PHeurMQLib(scip_, &graph_);
    // SCIP_CALL_EXC(SCIPincludeObjHeur(scip_, pHeurMqLib, TRUE))

    auto pHeurMst = new PHeurMst(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjHeur(scip_, pHeurMst, TRUE))

    auto pHeurKl = new PHeurKl(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjHeur(scip_, pHeurKl, TRUE))

    auto pHeurM = new PHeurMQLib(scip_, &graph_);
    SCIP_CALL_EXC(SCIPincludeObjHeur(scip_, pHeurM, TRUE))

    //auto pHeurB = new PHeurBurer(scip_, &graph_, 0);
    //SCIP_CALL_EXC(SCIPincludeObjHeur(scip_, pHeurB, TRUE))

    // Add branching heuristic
    auto brachingDegreeBased = new BranchruleDegree(scip_);
    SCIP_CALL_EXC(SCIPincludeObjBranchrule(scip_, brachingDegreeBased, TRUE))

    auto *disp = new ExtendedDisp(scip_);
    SCIP_CALL_EXC(SCIPincludeObjDisp(scip_, disp, true));

    addWarmStartSolutionsToScip();
    setScipParams();

    // Run the solver
    SCIP_CALL_EXC(SCIPsolve(scip_));

    if (!statsScip.empty()) {
        scipWriteStatistics(scip_, statsScip);
    }

    if (!history.empty() && eventHandler != nullptr) {
        auto eventHandlerLog = std::ofstream(history);
        eventHandler->to_json(eventHandlerLog);
        eventHandlerLog.close();
    }

    hasRun_ = true;
    std::cout << "Solver finished." << std::endl;

    SCIPprintStatistics(scip_, nullptr);

    return getBestSolution();
}

// set default params
void BasicOcwSolver::setScipParams() {
    // display
    SCIPsetIntParam(scip_, "display/freq", 1);
    SCIPsetIntParam(scip_, "display/currows/active", 2);            // 2 is on
    SCIPsetIntParam(scip_, "display/curdualbound/active", 2);       // 2 is on
    SCIPsetIntParam(scip_, "display/conflicts/active", 0);          // 2 is on
    SCIPsetIntParam(scip_, "display/vars/active", 0);               // 2 is on
    SCIPsetIntParam(scip_, "display/conss/active", 0);              // 2 is on
    SCIPsetIntParam(scip_, "display/curconss/active", 0);           // 2 is on
    SCIPsetIntParam(scip_, "display/curcols/active", 0);            // 2 is on
    SCIPsetIntParam(scip_, "display/nexternbranchcands/active", 0); // 2 is on
    SCIPsetIntParam(scip_, "display/depth/active", 2);              // 2 is on
    SCIPsetIntParam(scip_, "display/nfrac/active", 2);              // 2 is on
    SCIPsetIntParam(scip_, "display/separounds/active", 2);         // 2 is on
    SCIPsetIntParam(scip_, "display/poolsize/active", 2);           // 2 is on
    SCIPsetIntParam(scip_, "display/width", 180);                   // default 143
    // threads and lp solving
    SCIPsetIntParam(scip_, "lp/threads", 1);
    SCIPsetCharParam(scip_, "lp/initalgorithm", 'd');
    SCIPsetCharParam(scip_, "lp/resolvealgorithm", 'd');
    SCIPsetIntParam(scip_, "parallel/maxnthreads", 1);
    // fine tuned cutting
    SCIPsetRealParam(scip_, "separating/minefficacyroot", 0.000001); // default: 0.0001
    // SCIPsetIntParam(scip_, "limits/maxsol", 100);
    //  SCIPsetCharParam(scip_, "separating/orthofunc", 'd'); // default: 0.0001
    //  SCIPsetCharParam(scip_, "separating/efficacynorm", 'd'); // default: 0.0001
    //  SCIPsetBoolParam(scip_, "constraints/indicator/forcerestart", true); // default: 0.0001
    // SCIPsetIntParam(scip_, "limits/autorestartnodes", 500);
    SCIPsetIntParam(scip_, "separating/cutagelimit", 1);
    SCIPsetBoolParam(scip_, "separating/filtercutpoolrel", true);
    // SCIPsetRealParam(scip_, "cutselection/hybrid/minorthoroot", 0.5);
    SCIPsetIntParam(scip_, "separating/maxstallroundsroot", 600);
    // SCIPsetIntParam(scip_, "cutselection/dynamic/priority", 9000);
    //  SCIPsetIntParam(scip_, "separating/maxcutsroot", 10000);
    SCIPsetIntParam(scip_, "separating/poolfreq", 5);
    // load param file last, to overwrite the just set default values
    if (!paramFile_.empty()) {
        std::cout << "Using param file " << paramFile_ << std::endl;
        auto readReturn = SCIPreadParams(scip_, paramFile_.c_str());
        if (readReturn != SCIP_OKAY) {
            throw std::runtime_error("Error reading parameter file: " + paramFile_
                                     + ".\nSCIP return code: " + std::to_string(readReturn));
        }
    } else {
        std::cout << "No param file provided. Purely using default parameters." << std::endl;
    }
    // explicitly override some parameters
    SCIPsetIntParam(scip_, "randomization/randomseedshift", seed_);
    SCIPsetRealParam(scip_, "limits/time", timelimit_);
    SCIPsetLongintParam(scip_, "limits/totalnodes", nodelimit_);
}

void BasicOcwSolver::addWarmStartSolution(const std::vector<bool> &solution) {
    warmStartSolutions_.push_back(solution);
}

void BasicOcwSolver::addWarmStartSolutionsEdgesOnlyToScip() {
    for (auto &solution : warmStartSolutions_) {
        auto *probData = SCIPgetObjProbData(scip_);
        auto probDataMc = dynamic_cast<ProbDataEdgesMc *>(probData);

        SCIP_SOL *newSol;
        SCIP_CALL_EXC(SCIPcreateSol(scip_, &newSol, nullptr));

        for (auto e : graph_.edgeRange()) {
            SCIP_Real edgeVarValue = solution[e.u] ^ solution[e.v];
            SCIP_CALL_EXC(SCIPsetSolVal(scip_, newSol, probDataMc->edgeToVar(e.u, e.v), edgeVarValue));
        }

        assert(!SCIPsolIsPartial(newSol));

        SCIP_Bool added;
        SCIP_CALL_EXC(SCIPaddSolFree(scip_, &newSol, &added));
        assert(added);
    }
}

void BasicOcwSolver::addWarmStartSolutionsToScip() {
    for (auto &solution : warmStartSolutions_) {
        auto *probData = SCIPgetObjProbData(scip_);
        auto probDataRootedMc = dynamic_cast<ProbDataRootedMc *>(probData);
        auto root = probDataRootedMc->getRoot();

        if (solution[root]) {
            solution.flip();
        }

        SCIP_SOL *newSol;
        SCIP_CALL_EXC(SCIPcreateSol(scip_, &newSol, nullptr));

        for (unsigned int u = 0; u < graph_.upperNodeIdBound(); u++) {
            if (u != root) {
                SCIP_VAR *var = probDataRootedMc->nodeToVar(u);
                SCIP_CALL_EXC(SCIPsetSolVal(scip_, newSol, var, solution[u]));
            }
        }

        for (auto e : graph_.edgeRange()) {
            SCIP_Real edgeVarValue = solution[e.u] ^ solution[e.v];
            SCIP_CALL_EXC(SCIPsetSolVal(scip_, newSol, probDataRootedMc->edgeToVar(e.u, e.v), edgeVarValue));
        }

        assert(!SCIPsolIsPartial(newSol));

        SCIP_Bool added;
        SCIP_CALL_EXC(SCIPaddSolFree(scip_, &newSol, &added));
        assert(added);
    }
}

std::vector<bool> BasicOcwSolver::getBestSolution() const {
    assert(hasRun_);
    if (!nodeVars_)
        return getBestSolutionEdgesOnly();
    std::vector<bool> bestSolution(graph_.upperNodeIdBound(), false);

    auto *probData = SCIPgetObjProbData(scip_);
    auto probDataRootedMc = dynamic_cast<ProbDataRootedMc *>(probData);

    node root = probDataRootedMc->getRoot();
    bestSolution[root] = false;

    auto sol = SCIPgetBestSol(scip_);
    for (auto u : graph_.nodeRange()) {
        if (u != root) {
            bestSolution[u] = !SCIPisZero(scip_, SCIPgetSolVal(scip_, sol, probDataRootedMc->nodeToVar(u)));
        }
    }

    return bestSolution;
}

std::vector<bool> BasicOcwSolver::getBestSolutionEdgesOnly() const {
    assert(hasRun_);
    auto sol = SCIPgetBestSol(scip_);
    std::vector<Partition> bestSolution(graph_.upperNodeIdBound(), kUNASSIGNED);

    auto *probData = SCIPgetObjProbData(scip_);
    auto probDataMc = dynamic_cast<ProbDataEdgesMc *>(probData);

    std::vector<node> stack;

    for (unsigned int j = 0; j < graph_.upperNodeIdBound(); j++) {
        if (bestSolution[j] == sms::Partition::kUNASSIGNED) {
            stack.push_back(j);
            bestSolution[j] = kZERO;

            while (!stack.empty()) {
                auto currentNode = stack.back();
                stack.pop_back();

                for (auto curNeighbor : (&graph_)->neighborRange(currentNode)) {
                    auto edgeVar = probDataMc->edgeToVar(currentNode, curNeighbor);
                    assert(isBinary(scip_, edgeVar, sol));
                    auto edgeVal = SCIPgetSolVal(scip_, sol, edgeVar);
                    uint8_t binaryEdgeVal = edgeVal < 0.5 ? 0 : 1;

                    if (bestSolution[curNeighbor] == sms::Partition::kUNASSIGNED) {
                        if (binaryEdgeVal ^ (bestSolution[j] == kONE)) {
                            bestSolution[curNeighbor] = kONE;
                        } else {
                            bestSolution[curNeighbor] = kZERO;
                        }
                        stack.push_back(curNeighbor);
                    }
                }
            }
        }
    }

    std::vector<bool> bestSolutionBool(graph_.upperNodeIdBound(), false);
    for (unsigned int u = 0; u < graph_.upperNodeIdBound(); u++) {
        assert(bestSolution[u] == kZERO || bestSolution[u] == kONE);
        bestSolutionBool[u] = bestSolution[u];
    }
    return bestSolutionBool;
}

double BasicOcwSolver::bestSolutionValue() const {
    assert(hasRun_);
    return SCIPgetSolOrigObj(scip_, SCIPgetBestSol(scip_));
}

bool BasicOcwSolver::optimalityProven() const {
    return SCIPgetStatus(scip_) == SCIP_STATUS_OPTIMAL;
}

nlohmann::ordered_json BasicOcwSolver::getStats() const {
    nlohmann::ordered_json stats = {{"number of vertices", graph_.numberOfNodes()},
                                    {"number of edges", graph_.numberOfEdges()}};
    auto otherStats = scipFinalStatsToJson(scip_);
    // merge stats and other stats
    stats.merge_patch(otherStats);
    return stats;
}

} // namespace sms
