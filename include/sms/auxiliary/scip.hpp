#ifndef SMS_SCIP_HPP
#define SMS_SCIP_HPP

#include <string>

#include "nlohmann/json.hpp"
#include <scip/scip.h>

nlohmann::ordered_json scipFinalStatsToJson(SCIP *scip);

// make sure scip does not scale the objective function and assumes all solution values to be integer
void setObjectiveIntegral(SCIP *scip);

bool scipVarIsFixedLocal(SCIP *scip, SCIP_Var *var);

SCIP_Real normalizedFractionality(SCIP_Real);

bool consArrayIsValid(int nconss, SCIP_CONS **conss);

bool isBinary(Scip *scip, SCIP_Var *var);

bool isBinary(SCIP *scip, SCIP_VAR *var, SCIP_Sol *sol);

void scipWriteStatistics(SCIP *scip, const std::string &outputFileName);

void scipWriteBestSolution(SCIP *scip, const std::string &outputFileName);

bool isFreeLocal(SCIP_VAR *var);

bool isFreeGlobal(SCIP_VAR *var);

int numberOfLpThreads(SCIP *scip);

void addConstraintNoCapture(SCIP *scip, SCIP_Real coeffs[], SCIP_VAR *vars[], size_t nvars, SCIP_Real lhs,
                            SCIP_Real rhs);

void addXorConstraintNoCapture(SCIP *scip, SCIP_VAR *vars[], size_t nvars, const std::string &name);

void addXorConstraintNoCapture(SCIP *scip, SCIP_VAR *vars[], size_t nvars);

std::string scipVersion();

std::string lpSolver();

#endif // SMS_SCIP_HPP
