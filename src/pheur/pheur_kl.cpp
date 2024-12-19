#include "sms/pheur/pheur_kl.hpp"

#include <iostream>

#include "sms/instance/mc_solution.hpp"
#include "sms/probdata/mc_rooted.hpp"

#define HEUR_NAME "kl"
#define HEUR_DESC "Kernighan-Lin Heuristic for MaxCut"
#define HEUR_DISPCHAR 'K'
#define HEUR_PRIORITY (1)
#define HEUR_FREQ 1
#define HEUR_FREQOFS 0
#define HEUR_MAXDEPTH (-1)
#define HEUR_TIMING SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP FALSE

namespace sms {

PHeurKl::PHeurKl(SCIP *scip, const NetworKit::Graph *g)
    : ObjHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS, HEUR_MAXDEPTH,
              HEUR_TIMING, HEUR_USESSUBSCIP),
      mcGraph_(g),
      heuristic_(mcGraph_) {}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_RETCODE PHeurKl::scip_free(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_RETCODE PHeurKl::scip_init(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_RETCODE PHeurKl::scip_exit(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_RETCODE PHeurKl::scip_initsol(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_RETCODE PHeurKl::scip_exitsol(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

SCIP_RETCODE
PHeurKl::scip_exec(SCIP *scip, SCIP_HEUR *heur, SCIP_HEURTIMING /*heurtiming*/, unsigned int /*nodeinfeasible*/,
                   SCIP_RESULT *result) {
    assert(scip != nullptr);
    assert(result != nullptr);
    assert(heur != nullptr);

    *result = SCIP_DIDNOTRUN;

    SCIP_SOL *currentBestSol = SCIPgetBestSol(scip);
    if (currentBestSol == nullptr) {
        *result = SCIP_DIDNOTFIND;
        return SCIP_OKAY;
    }

    SCIP_SOL **solutions = SCIPgetSols(scip);
    auto numSols = SCIPgetNSols(scip);
    for (int i = 0; i < numSols; ++i) {
        SCIP_RESULT currentResult = SCIP_DIDNOTFIND;
        auto *nextSol = solutions[i];
        auto idNextSol = SCIPsolGetIndex(nextSol);
        if (idNextSol > lastSolutionIndex_ && (heur != SCIPgetSolHeur(scip, nextSol))) {
            improveSolution(scip, heur, nextSol, &currentResult);
            if (currentResult == SCIP_FOUNDSOL)
                *result = SCIP_FOUNDSOL;
        }
        lastSolutionIndex_ = std::max(idNextSol, lastSolutionIndex_);
    }

    return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
scip::ObjCloneable *PHeurKl::clone(SCIP *scip) const {
    return new PHeurKl(scip, mcGraph_);
}

SCIP_RETCODE PHeurKl::improveSolution(SCIP *scip, SCIP_HEUR *heur, SCIP_SOL *sol, SCIP_RESULT *result) {
    // update weights of heuristic
    auto *probData = SCIPgetObjProbData(scip);
    auto pProbDataRootedMc = dynamic_cast<ProbDataRootedMc *>(probData);
    assert(pProbDataRootedMc != nullptr);
    std::vector<bool> initSol(mcGraph_->upperNodeIdBound(), false);
    std::vector<NetworKit::node> fixedVertices;
    for (auto u : mcGraph_->nodeRange()) {
        if (u != pProbDataRootedMc->getRoot()) {
            SCIP_VAR *vertexVar = pProbDataRootedMc->nodeToVar(u);
            initSol[u] = SCIPgetSolVal(scip, sol, vertexVar) > 0.5;
            if (!isFreeGlobal(vertexVar))
                fixedVertices.push_back(u);
        }
    }
    initSol[pProbDataRootedMc->getRoot()] = false;
    heuristic_.phase1Optimization(initSol, fixedVertices);

    // add solution to scip
    auto solution = heuristic_.getPrimalSolution();

    SCIP_SOL *newSol;
    SCIP_CALL(SCIPcreateSol(scip, &newSol, heur));
    // root needs to be in partition zero
    if (solution[pProbDataRootedMc->getRoot()])
        solution.flip();
    for (unsigned int u = 0; u < mcGraph_->upperNodeIdBound(); u++) {
        if (pProbDataRootedMc->getRoot() != u)
            SCIP_CALL(SCIPsetSolVal(scip, newSol, pProbDataRootedMc->nodeToVar(u), solution[u]));
    }
    for (auto e : mcGraph_->edgeRange()) {
        SCIP_Real edgeVarValue = solution[e.u] ^ solution[e.v];
        SCIP_CALL(SCIPsetSolVal(scip, newSol, pProbDataRootedMc->edgeToVar(e.u, e.v), edgeVarValue));
    }
    assert(!SCIPsolIsPartial(newSol));
    assert(!SCIPsolIsOriginal(newSol));

    SCIP_Bool success;
    // try heuristic solution AND ADD IT to SCIP solution list (SCIP function name is misleading)
    SCIP_CALL(SCIPtrySolFree(scip, &newSol, TRUE, FALSE, FALSE, FALSE, TRUE, &success));
    // assert(success);

    if (success) {
        *result = SCIP_FOUNDSOL;
    }
    return SCIP_OKAY;
}

} // namespace sms
