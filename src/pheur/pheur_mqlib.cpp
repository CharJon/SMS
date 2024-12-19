#include "sms/pheur/pheur_mqlib.hpp"

#include <iostream>

#include "sms/instance/mc_solution.hpp"
#include "sms/probdata/mc_rooted.hpp"

#define HEUR_NAME "bur02"
#define HEUR_DESC "MQlib BURER2002 for MaxCut"
#define HEUR_DISPCHAR 'M'
#define HEUR_PRIORITY (1000000)
#define HEUR_FREQ 1
#define HEUR_FREQOFS 0
#define HEUR_MAXDEPTH 0
#define HEUR_TIMING SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP FALSE

#define DEFAULT_MAXSECONDS 2.0

namespace sms {

PHeurMQLib::PHeurMQLib(SCIP *scip, const Graph *g)
    : ObjHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS, HEUR_MAXDEPTH,
              HEUR_TIMING, FALSE),
      mcGraph_(g),
      heurRuntime_(DEFAULT_MAXSECONDS) {
    SCIP_CALL_EXC_NO_THROW(SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxseconds",
                                            "maximal number of seconds to run the heuristic", &heurRuntime_, FALSE,
                                            DEFAULT_MAXSECONDS, 0.0, SCIP_REAL_MAX, nullptr, nullptr));
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_RETCODE PHeurMQLib::scip_free(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_RETCODE PHeurMQLib::scip_init(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_RETCODE PHeurMQLib::scip_exit(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_RETCODE PHeurMQLib::scip_initsol(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_RETCODE PHeurMQLib::scip_exitsol(SCIP * /*scip*/, SCIP_HEUR * /*heur*/) {
    return SCIP_OKAY;
}

SCIP_RETCODE
PHeurMQLib::scip_exec(SCIP *scip, SCIP_HEUR *heur, SCIP_HEURTIMING /*heurtiming*/, unsigned int /*nodeinfeasible*/,
                      SCIP_RESULT *result) {
    assert(scip != nullptr);
    assert(result != nullptr);
    assert(heur != nullptr);

    *result = SCIP_DIDNOTFIND;

    auto currentCutoffBound = SCIPgetDualbound(scip);
    currentCutoffBound =
        std::floor(currentCutoffBound); // this is optimistic, the objective might not always have integer values

    auto *probData = SCIPgetObjProbData(scip);
    auto probDataRootedMc = dynamic_cast<ProbDataRootedMc *>(probData);

    HeuristicSolver heurSolver(*mcGraph_);
    heurSolver.setTimelimit(std::chrono::duration<double>(heurRuntime_));
    heurSolver.setSolutionValueLimit(currentCutoffBound);
    heurSolver.run();
    auto smsMcSol = heurSolver.getBestSolution();

    if (smsMcSol[probDataRootedMc->getRoot()])
        smsMcSol.flip();

    SCIP_SOL *newSol;
    SCIP_CALL(SCIPcreateSol(scip, &newSol, heur));

    for (unsigned int u = 0; u < mcGraph_->upperNodeIdBound(); u++) {
        if (u != probDataRootedMc->getRoot()) {
            SCIP_CALL(SCIPsetSolVal(scip, newSol, probDataRootedMc->nodeToVar(u), smsMcSol[u]));
        }
    }

    for (auto e : mcGraph_->edgeRange()) {
        SCIP_Real edgeVarValue = smsMcSol[e.u] ^ smsMcSol[e.v];
        SCIP_CALL(SCIPsetSolVal(scip, newSol, probDataRootedMc->edgeToVar(e.u, e.v), edgeVarValue));
    }

    assert(!SCIPsolIsPartial(newSol));
    assert(!SCIPsolIsOriginal(newSol));

    SCIP_Bool success;
    // try heuristic solution AND ADD IT to SCIP solution list (SCIP function name is misleading)
    SCIP_CALL(SCIPtrySolFree(scip, &newSol, TRUE, FALSE, FALSE, FALSE, TRUE, &success));

    if (success) {
        *result = SCIP_FOUNDSOL;
    }

    return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
scip::ObjCloneable *PHeurMQLib::clone(SCIP *scip) const {
    return new PHeurMQLib(scip, mcGraph_);
}

} // namespace sms
