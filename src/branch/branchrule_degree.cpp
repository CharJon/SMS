#include "sms/branch/branchrule_degree.hpp"

#include "sms/probdata/mc_rooted.hpp"

#define BRANCHRULE_NAME "degree"
#define BRANCHRULE_DESC "branches on the vertex-variable having the highest (dynamic) degree"
#define BRANCHRULE_PRIORITY 30000   // highest default is 10000
#define BRANCHRULE_MAXDEPTH (-1)    // -1 means no maximum
#define BRANCHRULE_MAXBOUNDDIST 1.0 // 1 means always
#define BRANCHRULE_DYNAMIC true     // default value
#define BRANCHRULE_TIEBREAK false   // default value

BranchruleDegree::BranchruleDegree(SCIP *scip)
    : scip::ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                          BRANCHRULE_MAXBOUNDDIST),
      dynamic_(BRANCHRULE_DYNAMIC),
      tieBreak_(BRANCHRULE_TIEBREAK) {
    SCIP_CALL_EXC_NO_THROW(SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/dynamic",
                                            "should degree be calculated based on free variables only", &dynamic_,
                                            FALSE, BRANCHRULE_DYNAMIC, nullptr, nullptr));
    SCIP_CALL_EXC_NO_THROW(SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/tiebreak", "should ties be broken",
                                            &tieBreak_, FALSE, BRANCHRULE_TIEBREAK, nullptr, nullptr));
}

/***
 * Calculates the (dynamic) degree of the corresponding vertex of a variable.
 * @param variable
 * @return the calculated degree
 */
int BranchruleDegree::degree(SCIP_VAR *variable) {
    auto *probData = SCIPgetObjProbData(scip_);
    auto pProbDataRootedMc = dynamic_cast<ProbDataRootedMc *>(probData);
    assert(pProbDataRootedMc != nullptr);

    auto root = pProbDataRootedMc->getRoot();
    auto u = pProbDataRootedMc->varToNode(variable);

    uint64_t degree = pProbDataRootedMc->getMcGraph()->degree(u);
    if (!dynamic_)
        return degree;
    for (auto neighbor : pProbDataRootedMc->getMcGraph()->neighborRange(u)) {
        if ((neighbor == root) || isBinary(scip_, pProbDataRootedMc->nodeToVar(neighbor))) {
            degree -= 1;
        }
    }

    return degree;
}

SCIP_RETCODE BranchruleDegree::scip_execlp(SCIP *scip, SCIP_BRANCHRULE * /*branchrule*/, unsigned int /*allowaddcons*/,
                                           SCIP_RESULT *result) {
    *result = SCIP_DIDNOTRUN;

    // get LP branching candidates
    int number_of_cands;
    SCIP_VAR **candidates;
    SCIP_Real *candidateValues;
    SCIP_CALL(SCIPgetLPBranchCands(scip, &candidates, &candidateValues, nullptr, nullptr, &number_of_cands, nullptr));

    uint64_t bestDegreeScore = 0;
    SCIP_Real bestPseudoCost = 0;
    SCIP_Real bestRootDiff = 0;

    int bestCandIdx = 0;

    for (int i = 0; i < number_of_cands; i++) {
        SCIP_VAR *currentVar = candidates[i];
        auto currentVal = candidateValues[i];
        uint64_t currentDegree = degree(currentVar);
        SCIP_Real currentPseudocost = SCIPgetVarPseudocostScore(scip, currentVar, currentVal);
        SCIP_Real currentRootSolVal = SCIPvarGetRootSol(currentVar);
        SCIP_Real currentRootDiff = REALABS(currentVal - currentRootSolVal);

        if ((currentDegree > bestDegreeScore)
            || (tieBreak_ && (currentDegree == bestDegreeScore) && SCIPisSumGT(scip, currentPseudocost, bestPseudoCost))
            || (tieBreak_ && (currentDegree == bestDegreeScore) && SCIPisSumEQ(scip, currentPseudocost, bestPseudoCost)
                && (currentRootDiff > bestRootDiff))) {
            bestDegreeScore = currentDegree;
            bestPseudoCost = currentPseudocost;
            bestRootDiff = currentRootDiff;
            bestCandIdx = i;
        }
    }

    SCIP_NODE *downchild;
    SCIP_NODE *upchild;
    SCIP_NODE *eqchild;

    SCIPbranchVar(scip_, candidates[bestCandIdx], &downchild, &eqchild, &upchild);
    assert(downchild != nullptr);
    assert(upchild != nullptr);

    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}
