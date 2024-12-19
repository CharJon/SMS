#include "scip/scip.h"

#include <fstream>
#include <string>

#include "nlohmann/json.hpp"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_xor.h"

#include "sms/auxiliary/scip.hpp"
#include "sms/auxiliary/scip_exception.hpp"

void addConstraintNoCapture(SCIP *scip, SCIP_Real coeffs[], SCIP_VAR *vars[], size_t nvars, SCIP_Real lhs,
                            SCIP_Real rhs) {
    SCIP_CONS *cons;
    SCIP_CALL_EXC(SCIPcreateConsBasicLinear(scip, &cons, "", nvars, vars, coeffs, lhs, rhs))
    SCIP_CALL_EXC(SCIPaddCons(scip, cons))
    SCIP_CALL_EXC(SCIPreleaseCons(scip, &cons))
}

void addLogicOrConstraintNoCapture(SCIP *scip, SCIP_VAR *vars[], size_t nvars) {
    SCIP_CONS *cons;
    SCIP_CALL_EXC(SCIPcreateConsBasicLogicor(scip, &cons, "", nvars, vars))
    SCIP_CALL_EXC(SCIPaddCons(scip, cons))
    SCIP_CALL_EXC(SCIPreleaseCons(scip, &cons))
}

void addXorConstraintNoCapture(SCIP *scip, SCIP_VAR *vars[], size_t nvars) {
    addXorConstraintNoCapture(scip, vars, nvars, "");
}

void addXorConstraintNoCapture(SCIP *scip, SCIP_VAR *vars[], size_t nvars, const std::string &name) {
    SCIP_CONS *cons;
    SCIP_CALL_EXC(SCIPcreateConsXor(scip, &cons, name.c_str(), false, nvars, vars, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                                    FALSE, FALSE, FALSE, FALSE))
    SCIP_CALL_EXC(SCIPaddCons(scip, cons))
    SCIP_CALL_EXC(SCIPreleaseCons(scip, &cons))
}

nlohmann::ordered_json scipFinalStatsToJson(SCIP *scip) {
    nlohmann::ordered_json j;
    j["mip solver"] = scipVersion();
    j["lp solver"] = lpSolver();
    j["cpu_solving_time"] = SCIPgetSolvingTime(scip);
    j["best_solution_value"] = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
    j["optimal_solution_found"] = SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL;
    j["timelimit_reached"] = SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT;
    j["bab_nodes"] = SCIPgetNTotalNodes(scip);
    j["final_bound_value"] = SCIPgetDualbound(scip);
    j["root_bound_value"] = SCIPgetDualboundRoot(scip);
    int seed = 0;
    SCIPgetIntParam(scip, "randomization/randomseedshift", &seed);
    j["scip_randomseedshift"] = seed;
    double timelimit = 0;
    SCIPgetRealParam(scip, "limits/time", &timelimit);
    j["scip_timelimit"] = timelimit;
    return j;
}

bool consArrayIsValid(int nconss, SCIP_CONS **conss) {
    bool isValid = true;
    for (int i = 0; i < nconss; ++i) {
        isValid &= (conss != nullptr);
        isValid &= (conss[i] != nullptr);
    }
    return isValid;
}

bool scipVarIsFixedLocal(SCIP *scip, SCIP_Var *var) {
    return SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
}

// non continuous variables only!
bool isFreeLocal(SCIP_VAR *var) {
    // assert(var->vartype != SCIP_VARTYPE_CONTINUOUS);
    auto diff = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
    return diff > 0.5;
}

// non continuous variables only!
bool isFreeGlobal(SCIP_VAR *var) {
    // assert(var->vartype != SCIP_VARTYPE_CONTINUOUS);
    auto diff = SCIPvarGetUbGlobal(var) - SCIPvarGetLbGlobal(var);
    return diff > 0.5;
}

void scipWriteStatistics(SCIP *scip, const std::string &outputFileName) {
    auto file = fopen(outputFileName.c_str(), "w");
    SCIPprintStatistics(scip, file);
    fclose(file);
}

void scipWriteBestSolution(SCIP *scip, const std::string &outputFileName) {
    auto file = fopen(outputFileName.c_str(), "w");
    SCIPprintBestSol(scip, file, FALSE);
    fclose(file);
}

bool isBinary(SCIP *scip, SCIP_VAR *var, SCIP_Sol *sol) {
    auto lpVal = SCIPgetSolVal(scip, sol, var);
    bool isZero = SCIPisEQ(scip, lpVal, 0);
    bool isOne = SCIPisEQ(scip, lpVal, 1);
    return isZero || isOne;
}

bool isBinary(SCIP *scip, SCIP_VAR *var) {
    auto lpVal = SCIPgetSolVal(scip, nullptr, var);
    bool isZero = SCIPisEQ(scip, lpVal, 0);
    bool isOne = SCIPisEQ(scip, lpVal, 1);
    return isZero || isOne;
}

int numberOfLpThreads(SCIP *scip) {
    int nThreads = -1;
    SCIP_LPI *lpi = nullptr;
    SCIP_CALL_EXC(SCIPgetLPI(scip, &lpi))
    if (lpi != nullptr) {
        SCIP_CALL_EXC(SCIPlpiGetIntpar(lpi, SCIP_LPPAR_THREADS, &nThreads))
    }
    return nThreads;
}

void setObjectiveIntegral(SCIP *scip) {
    SCIPsetBoolParam(scip, "misc/scaleobj", FALSE);
    SCIPsetObjIntegral(scip);
}

std::string lpSolver() {
    return SCIPlpiGetSolverName();
}

std::string scipVersion() {
    return "SCIP " + std::to_string(SCIPmajorVersion()) + "." + std::to_string(SCIPminorVersion()) + "."
           + std::to_string(SCIPtechVersion());
}
