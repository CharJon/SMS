#include "sms/auxiliary/scip.hpp"
#include "sms/probdata/prob_data_mc.hpp"

namespace sms {

void buildGeneralModel(Scip *scip, const MaxCut &mc, bool allBinary) {
    assert(scip != nullptr);

    auto *probdata = new ProbDataMc(mc.getGraph());
    SCIP_CALL_EXC(SCIPcreateObjProb(scip, "mc default model", probdata, true))

    SCIP_CALL_EXC(SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE))

    // Partition vars for each node
    for (auto u : mc.getGraph().nodeRange()) {
        SCIP_VAR *var;
        std::stringstream varname;
        varname << "z_" << u;
        SCIP_CALL_EXC(SCIPcreateVarBasic(scip, &var, varname.str().c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY))
        SCIP_CALL_EXC(SCIPaddVar(scip, var))
        probdata->addVarToNode(u, var, scip);
        assert(SCIPvarGetNUses(var) == 3);
        SCIP_CALL_EXC(SCIPreleaseVar(scip, &var))
    }

    // Cut vars for each edge
    for (auto e : mc.getGraph().edgeWeightRange()) {
        assert(e.u != e.v);
        auto u = std::min(e.u, e.v);
        auto v = std::max(e.u, e.v);
        auto weight = e.weight;

        SCIP_VAR *var;
        std::stringstream varname;
        varname << "y_" << u << "_" << v;
        SCIP_CALL_EXC(SCIPcreateVarBasic(scip, &var, varname.str().c_str(), 0, 1, weight,
                                         allBinary ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_IMPLINT))
        SCIP_CALL_EXC(SCIPaddVar(scip, var))
        probdata->addVarToEdge(NetworKit::Edge{u, v}, var, scip);
        assert(SCIPvarGetNUses(var) == 3);
        SCIP_CALL_EXC(SCIPreleaseVar(scip, &var))
    }

    // Constraints
    for (auto e : mc.getGraph().edgeRange()) {
        SCIP_VAR *vars[3] = {probdata->nodeToVar(e.u), probdata->edgeToVar(e), probdata->nodeToVar(e.v)};
        addXorConstraintNoCapture(scip, vars, 3);
    }
}

ProbDataMc *getProbDataMc(SCIP *scip) {
    assert(scip != nullptr);
    scip::ObjProbData *probdata = SCIPgetObjProbData(scip);
    assert(probdata != nullptr);
    auto *mcProbData = dynamic_cast<ProbDataMc *>(probdata);
    assert(mcProbData != nullptr);
    return mcProbData;
}

} // namespace sms