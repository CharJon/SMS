#include "sms/conshdlr/conshdlr_cliques.hpp"

#include <cassert>
#include <iostream>

#include "objscip/objscip.h"
#include "scip/cons_linear.h"

#include "networkit/clique/MaximalCliques.hpp"
#include "networkit/graph/Graph.hpp"
#include <networkit/graph/GraphTools.hpp>

#include "sms/auxiliary/scip.hpp"
#include "sms/conshdlr/maximal_cliques.hpp"
#include "sms/graph/graphs_statistics.hpp"
#include "sms/probdata/mc_rooted.hpp"

/* constraint handler properties */
#define CONSHDLR_NAME "cliquecut"
#define CONSHDLR_DESC "clique constraints"
#define CONSHDLR_SEPAPRIORITY (+700000)  /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY (-600000)  /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY (-600000) /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ (13) /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ (1)  /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ (100)
/**< frequency for using all instead of only the useful constraints in separation,
 *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS                                                                                          \
    (-1) /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS FALSE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING SCIP_PRESOLTIMING_ALWAYS

#define CONSHDLR_SEPA_CUTPOOL TRUE   /**< should separated constraints be added to the cutpool? */
#define CONSHDLR_SEPA_LIMIT 400ul    /**< How many cliques shall be considered for separation? */
#define CONSHDLR_SUBCLIQUESIZE_UB 24 /**< Upper bound on subclique size */

namespace sms {

ConshdlrCliques::ConshdlrCliques(Scip *scip, const Graph &sepGraph)
    : scip::ObjConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC, CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY,
                        CONSHDLR_CHECKPRIORITY, CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ,
                        CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
                        CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING),
      sepGraph_(sepGraph),
      cliqueSeparator_(sepGraph),
      sepaCutpool_(CONSHDLR_SEPA_CUTPOOL),
      subCliqueSizeUpperBound_(CONSHDLR_SUBCLIQUESIZE_UB) {
    SCIP_CALL_EXC_NO_THROW(SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/sepacutpool",
                                            "separate cuts into cutpool", &sepaCutpool_, FALSE, CONSHDLR_SEPA_CUTPOOL,
                                            nullptr, nullptr));
    SCIP_CALL_EXC_NO_THROW(SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxsubliquesize",
                                           "upper bound on subclique size", &subCliqueSizeUpperBound_, FALSE,
                                           CONSHDLR_SUBCLIQUESIZE_UB, 3, std::numeric_limits<int>::max(), nullptr,
                                           nullptr));
}

SCIP_RETCODE
ConshdlrCliques::scip_check(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS ** /*conss*/, int /*nconss*/,
                            SCIP_SOL * /*sol*/, SCIP_Bool /*checkintegrality*/, SCIP_Bool /*checklprows*/,
                            SCIP_Bool /*printreason*/, SCIP_Bool /*completely*/, SCIP_RESULT *result) {
    assert(result != nullptr);
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_RETCODE
ConshdlrCliques::scip_sepasol(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS ** /*conss*/, int /*nconss*/,
                              int /*nusefulconss*/, SCIP_SOL * /*sol*/, SCIP_RESULT * /*result*/) {
    throw std::runtime_error("Not implemented: scip_sepasol.");
}

SCIP_RETCODE ConshdlrCliques::scip_sepalp(SCIP *scip, SCIP_CONSHDLR *conshdlr, SCIP_CONS **conss, int nconss,
                                          int nusefulconss, SCIP_RESULT *result) {
    SCIP_CALL(sepaCliquesLp(scip, conshdlr, conss, nconss, nusefulconss, result));
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCliques::sepaOneCliqueLp(SCIP *scip, SCIP_CONSHDLR *conshdlr, SCIP_RESULT *result,
                                              const std::vector<NetworKit::node> &left,
                                              const std::vector<NetworKit::node> &right) const {
    auto *probData = SCIPgetObjProbData(scip);
    auto probDataSunMc = dynamic_cast<ProbDataEdgesMc *>(probData);
    assert(probDataSunMc != nullptr);

    auto leftSizeD = static_cast<double>(left.size());
    auto rightSizeD = static_cast<double>(right.size());
    auto cliqueSizeD = leftSizeD + rightSizeD;
    double rhs = floor(cliqueSizeD / 2.) * ceil(cliqueSizeD / 2.) - leftSizeD * rightSizeD;
    SCIP_Real lhs = -SCIPinfinity(scip);

    SCIP_ROW *row;
    SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "c", lhs, rhs, FALSE, FALSE, TRUE));

    SCIP_CALL(SCIPcacheRowExtensions(scip, row));
    for (auto u : left) {
        for (auto v : right) {
            SCIP_CALL(SCIPaddVarToRow(scip, row, probDataSunMc->edgeToVar(u, v), -1));
        }
    }

    for (unsigned int i = 0; i < left.size(); i++) {
        for (unsigned int j = i + 1; j < left.size(); j++) {
            SCIP_CALL(SCIPaddVarToRow(scip, row, probDataSunMc->edgeToVar(left[i], left[j]), 1));
        }
    }

    for (unsigned int i = 0; i < right.size(); i++) {
        for (unsigned int j = i + 1; j < right.size(); j++) {
            SCIP_CALL(SCIPaddVarToRow(scip, row, probDataSunMc->edgeToVar(right[i], right[j]), 1));
        }
    }
    SCIP_CALL(SCIPflushRowExtensions(scip, row));

    // add row
    if (SCIPisCutEfficacious(scip, nullptr, row)) {
        // auto efficacy = SCIPgetCutEfficacy(scip, nullptr, row);
        // std::cout << "Efficacy: " << efficacy << "\n";
        *result = SCIP_SEPARATED;
        if (sepaCutpool_) {
            SCIP_CALL(SCIPaddPoolCut(scip, row));
        } else {
            SCIP_Bool infeasible;
            SCIP_CALL(SCIPaddRow(scip, row, false, &infeasible));
            if (infeasible) {
                *result = SCIP_CUTOFF;
            }
        }
    }
    SCIP_CALL(SCIPreleaseRow(scip, &row));

    return SCIP_OKAY;
}

SCIP_RETCODE
ConshdlrCliques::sepaCliquesLp(SCIP *scip,              /**< SCIP data structure */
                               SCIP_CONSHDLR *conshdlr, /**< the constraint handler itself */
                               SCIP_CONS ** /*conss*/,  /**< array of constraints to process */
                               int /*nconss*/,          /**< number of constraints to process */
                               int /*nusefulconss*/,    /**< number of useful (non-obsolete) constraints to process */
                               SCIP_RESULT *result) {
    assert(scip != nullptr);
    *result = SCIP_DIDNOTFIND;

    auto *probData = SCIPgetObjProbData(scip);
    auto probDataSunMc = dynamic_cast<ProbDataRootedMc *>(probData);
    assert(probDataSunMc != nullptr);

    // Update lp values
    for (auto e : sepGraph_.edgeRange()) {
        auto var = probDataSunMc->edgeToVar(e.u, e.v);
        auto val = SCIPgetSolVal(scip, nullptr, var);
        auto fixedVal = std::min(std::max(0., val), 1.);
        cliqueSeparator_.updateWeights(e.u, e.v, fixedVal);
    }

    // loop over cliques and test for violated constraints
    bool cutoff = false;
    int numAddedThisRound = 0;
    for (unsigned int i = 0; (!cutoff) && (i < cliqueSeparator_.numCliques()); i++) {
        if (fails_[i] > 5) {
            continue;
        }
        // only process all cliques, if the first rounds did not find any cuts
        if ((i >= CONSHDLR_SEPA_LIMIT) && (numAddedThisRound > 0)) {
            break;
        }

        std::vector<NetworKit::node> currentClique;
        for (auto u : cliqueSeparator_.getClique(i)) {
            if (u == probDataSunMc->getRoot() || (!scipVarIsFixedLocal(scip, probDataSunMc->nodeToVar(u))))
                currentClique.push_back(u);
        }
        if (currentClique.size() <= 3)
            continue;
        auto violations = cliqueSeparator_.checkViolation(currentClique, 5, subCliqueSizeUpperBound_, true);
        bool addedAtLeastOne = false;
        for (const auto &res : violations) {
            SCIP_RESULT resSepa = SCIP_DIDNOTFIND;
            if ((res.first.size() + res.second.size()) != 0) {
                sepaOneCliqueLp(scip, conshdlr, &resSepa, res.first, res.second);
            }
            if (resSepa == SCIP_CUTOFF) {
                *result = SCIP_CUTOFF;
                cutoff = true;
                break;
            } else if (resSepa == SCIP_SEPARATED) {
                *result = SCIP_SEPARATED;
                numAddedThisRound++;
                addedAtLeastOne = true;
            }
        }
        if (not addedAtLeastOne) {
            fails_[i] += 1;
        }
    }

    return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_RETCODE ConshdlrCliques::scip_enfolp(SCIP *scip, SCIP_CONSHDLR *conshdlr, SCIP_CONS **conss, int nconss,
                                          int nusefulconss, unsigned int /*solinfeasible*/, SCIP_RESULT *result) {
    assert(result != nullptr);
    for (int i = 0; i < nconss; ++i) {
        assert(conss != nullptr);
        assert(conss[i] != nullptr);
    }

    *result = SCIP_INFEASIBLE; // Infeasible does not really mean infeasible here. See SCIP docu
    SCIP_RESULT resultSepa;
    SCIP_CALL(sepaCliquesLp(scip, conshdlr, conss, nconss, nusefulconss, &resultSepa));
    if (resultSepa == SCIP_SEPARATED || resultSepa == SCIP_NEWROUND)
        *result = SCIP_SEPARATED;
    else {
        *result = SCIP_FEASIBLE;
    }

    return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_RETCODE
ConshdlrCliques::scip_enfops(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS ** /*conss*/, int /*nconss*/,
                             int /*nusefulconss*/, SCIP_Bool /*solinfeasible*/, SCIP_Bool /*objinfeasible*/,
                             SCIP_RESULT *result) {
    *result = SCIP_FEASIBLE;

    throw std::runtime_error("Not implemented: scip_enfops");
}

SCIP_RETCODE
ConshdlrCliques::scip_lock(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS * /*cons*/,
                           SCIP_LOCKTYPE /*locktype*/, int /*nlockspos*/, int /*nlocksneg*/) {
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCliques::scip_trans(SCIP *scip, SCIP_CONSHDLR *conshdlr, SCIP_CONS *sourcecons,
                                         SCIP_CONS **targetcons) {

    // SCIP_CONSDATA *sourcedata;
    SCIP_CONSDATA *targetdata = nullptr;

    // create target constraint
    SCIP_CALL(SCIPcreateCons(
        scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata, SCIPconsIsInitial(sourcecons),
        SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
        SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
        SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)));

    return SCIP_OKAY;
}

scip::ObjProbCloneable *ConshdlrCliques::clone(SCIP *scip, SCIP_Bool *valid) const {
    assert(valid != nullptr);
    *valid = true;
    return new ConshdlrCliques(scip, sepGraph_);
}

SCIP_RETCODE
ConshdlrCliques::scip_initsol(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS ** /*conss*/, int /*nconss*/) {
    // collect all cliques of odd size with at least 5 nodes

    const int enumerationLimit = 10 * CONSHDLR_SEPA_LIMIT;
    int numCliquesEnumerated = 0;
    // stop if either enough cliques of size >= 5 were found or the enumeration limit is reached
    auto callback = [this, &numCliquesEnumerated](const std::vector<NetworKit::node> &clique) {
        numCliquesEnumerated += 1;
        if (clique.size() >= 5) {
            cliqueSeparator_.addClique(clique);
        }
        return numCliquesEnumerated >= enumerationLimit;
    };

    auto mc = MaximalCliques(sepGraph_, callback);
    mc.run();

    cliqueSeparator_.sortCliquesBySize();
    fails_ = std::vector<int>(cliqueSeparator_.numCliques(), 0);
    std::cout << "Cliques: " << cliqueSeparator_.numCliques() << std::endl;
    std::cout << "Large: " << cliqueSeparator_.numLargeCliques() << std::endl;

    return SCIP_OKAY;
}

} // namespace sms
