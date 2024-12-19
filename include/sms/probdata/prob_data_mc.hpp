#ifndef SMS_MC_EDGES_HPP
#define SMS_MC_EDGES_HPP

#include "objscip/objscip.h"
#include "scip/scipdefplugins.h"

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/scip_exception.hpp"
#include "sms/instance/maxcut.hpp"

namespace sms {

/*
 * Main purpose of this class:
 * Store the graph and mappings in both directions, for
 * Nodes <-> node variables
 * Edges <-> edge variables
 */
class ProbDataMc : public scip::ObjProbData {
public:
    explicit ProbDataMc(const NetworKit::Graph &g) : graph_(g), nodeIdToVar_(), edgeIdToVar_() {
        graph_.indexEdges();
        nodeIdToVar_.resize(graph_.upperNodeIdBound(), nullptr);
        edgeIdToVar_.resize(graph_.upperEdgeIdBound(), nullptr);
    }

    ~ProbDataMc() override = default; // to not forget override, if we actually need destructor

    void addVarToNode(NetworKit::node u, SCIP_VAR *var, SCIP *scip) {
        assert(graph_.hasNode(u));
        SCIP_CALL_EXC(SCIPcaptureVar(scip, var))
        nodeIdToVar_[u] = var;
        int varIndex = SCIPvarGetIndex(var);
        varIdToNode_[varIndex] = u;
    }

    void addVarToEdge(NetworKit::Edge e, SCIP_VAR *var, SCIP *scip) {
        assert(graph_.hasEdge(e.u, e.v));
        SCIP_CALL_EXC(SCIPcaptureVar(scip, var))
        auto edgeId = graph_.edgeId(e.u, e.v);
        edgeIdToVar_[edgeId] = var;
        int varIndex = SCIPvarGetIndex(var);
        varIdToEdge_[varIndex] = e;
    }

    void addEdge(SCIP *scip, NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        // add edge to graph
        graph_.addEdge(u, v, weight);
        // add variable and store it
        SCIP_Var *var = nullptr;
        std::stringstream varname;
        varname << "y_" << u << "_" << v;
        SCIP_CALL_EXC(SCIPcreateVarBasic(scip, &var, varname.str().c_str(), 0, 1, weight,
                                         SCIP_VARTYPE_BINARY)) // var gets captured here
        SCIP_CALL_EXC(SCIPaddVar(scip, var))
        NetworKit::edgeid edgeId = graph_.edgeId(u, v);
        if (edgeId >= edgeIdToVar_.size()) {
            edgeIdToVar_.resize(2 * edgeIdToVar_.size(), nullptr);
        }
        edgeIdToVar_[edgeId] = var;
        int varIndex = SCIPvarGetIndex(var);
        varIdToEdge_[varIndex] = {u, v};
        // ToDo (JC): What about constraints?!
    }

    SCIP_VAR *edgeToVar(NetworKit::Edge e) const { return edgeIdToVar(graph_.edgeId(e.u, e.v)); }

    SCIP_VAR *edgeIdToVar(NetworKit::index edgeId) const { return edgeIdToVar_[edgeId]; }

    SCIP_VAR *nodeToVar(NetworKit::node u) const { return nodeIdToVar_[u]; }

    // SCIP Methods

    /** destructor of user problem data to free original user data (called when original problem is freed)
     *
     *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
     *  because all the work to delete the user problem data can be done in the destructor of the user problem
     *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
     *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
     *  longer needed.
     */
    SCIP_RETCODE scip_delorig(SCIP *scip) override {
        freeCapturedVars(scip);
        return SCIP_OKAY;
    }

    /** creates user data of transformed problem by transforming the original user problem data
     *  (called after problem was transformed)
     *
     *  The user has two possibilities to implement this method:
     *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
     *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
     *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
     *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
     *      solving process is terminated.
     *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
     *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
     *      destructor of the object if the transformed problem data is no longer needed.
     */
    SCIP_RETCODE
    scip_trans(SCIP *scip,                /**< SCIP data structure */
               ObjProbData **objprobdata, /**< pointer to store the transformed problem data object */
               SCIP_Bool *deleteobject    /**< pointer to store whether SCIP should delete the object after solving */
               ) override {
        assert(objprobdata != nullptr);
        assert(deleteobject != nullptr);

        auto *propData = new ProbDataMc(graph_);
        *objprobdata = propData;
        *deleteobject = true;

        // set all variable mappings to the transformed vars
        for (auto u : graph_.nodeRange()) {
            SCIP_VAR *transvar = nullptr;
            SCIP_CALL_EXC_NO_THROW(SCIPgetTransformedVar(scip, nodeIdToVar_[u], &transvar))
            propData->addVarToNode(u, transvar, scip);
        }
        for (auto e : graph_.edgeRange()) {
            SCIP_VAR *transvar = nullptr;
            SCIP_CALL_EXC_NO_THROW(SCIPgetTransformedVar(scip, edgeIdToVar_[graph_.edgeId(e.u, e.v)], &transvar))
            propData->addVarToEdge(e, transvar, scip);
        }

        return SCIP_OKAY;
    }

    /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
     *
     *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
     *  because all the work to delete the user problem data can be done in the destructor of the user problem
     *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
     *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
     *  longer needed.
     */
    SCIP_RETCODE scip_deltrans(SCIP *scip /**< SCIP data structure */
                               ) override {
        freeCapturedVars(scip);
        return SCIP_OKAY;
    }

    /** solving process initialization method of transformed data (called before the branch and bound process begins)
     *
     *  This method is called before the branch and bound process begins and can be used to initialize user problem
     *  data that depends for example on the number of active problem variables, because these are now fixed.
     */
    SCIP_RETCODE scip_initsol(SCIP * /*scip*/ /**< SCIP data structure */
                              ) override {
        // ToDo (JC): Compactify graph here!
        return SCIP_OKAY;
    }

    /** solving process deinitialization method of transformed data (called before the branch and bound data is freed)
     *
     *  This method is called before the branch and bound data is freed and should be used to free all data that
     *  was allocated in the solving process initialization method. The user has to make sure, that all LP rows
     *  associated to the transformed user problem data are released.
     */
    SCIP_RETCODE scip_exitsol(SCIP * /*scip*/,      /**< SCIP data structure */
                              SCIP_Bool /*restart*/ /**< was this exit solve call triggered by a restart? */
                              ) override {
        return SCIP_OKAY;
    }

    /** copies user data of source SCIP for the target SCIP
     *
     *  This method should copy the problem data of the source SCIP and create a target problem data for (target)
     *  SCIP. Implementing this callback is optional. If the copying process was successful the target SCIP gets this
     *  problem data assigned. In case the result pointer is set to SCIP_DIDNOTRUN the target SCIP will have no problem
     * data at all.
     *
     *  The variable map and the constraint map can be used via the function SCIPgetVarCopy() and SCIPgetConsCopy(),
     *  respectively, to get for certain variables and constraints of the source SCIP the counter parts in the target
     *  SCIP. You should be very carefully in using these two methods since they could lead to an infinite loop due to
     *  recursion.
     *
     *  possible return values for *result:
     *  - SCIP_DIDNOTRUN  : the copying process was not performed
     *  - SCIP_SUCCESS    : the copying process was successfully performed
     */
    SCIP_RETCODE scip_copy(SCIP * /*scip*/,            /**< SCIP data structure */
                           SCIP * /*sourcescip*/,      /**< source SCIP main data structure */
                           SCIP_HASHMAP * /*varmap*/,  /**< a hashmap which stores the mapping of source
                                                        * variables to  corresponding  target variables */
                           SCIP_HASHMAP * /*consmap*/, /**< a hashmap which stores the mapping of source
                                                        * contraints to corresponding target constraints */
                           ObjProbData **objprobdata,  /**< pointer to store the copied problem data object */
                           SCIP_Bool /*global*/,       /**< create a global or a local copy? */
                           SCIP_RESULT *result         /**< pointer to store the result of the call */
                           ) override {
        (*objprobdata) = nullptr;
        (*result) = SCIP_DIDNOTRUN;
        throw std::runtime_error("Not implemented");
    }

    const NetworKit::Graph &getGraph() const { return graph_; }

private:
    NetworKit::Graph graph_;

    std::vector<SCIP_VAR *> nodeIdToVar_;
    std::unordered_map<int, NetworKit::node> varIdToNode_;

    std::vector<SCIP_VAR *> edgeIdToVar_;
    std::unordered_map<int, NetworKit::Edge> varIdToEdge_;

    void freeCapturedVars(SCIP *scip) {
        for (auto var : nodeIdToVar_) {
            // auto numUses = SCIPvarGetNUses(var);
            if (var != nullptr) {
                SCIP_CALL_EXC(SCIPreleaseVar(scip, &var))
            }
        }

        for (auto var : edgeIdToVar_) {
            // auto numUses = SCIPvarGetNUses(var);
            if (var != nullptr) {
                SCIP_CALL_EXC(SCIPreleaseVar(scip, &var))
            }
        }
    }

    void aggregateToEqualPartition(SCIP *scip, NetworKit::node nodeToRemove, NetworKit::node nodeToKeep,
                                   SCIP_Var *varToRemove, SCIP_Var *varToKeep) const {
        SCIP_Bool infeasible, redundant, aggregated, fixed;
        SCIPaggregateVars(scip, varToRemove, varToKeep, 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated);
        assert(SCIPvarGetStatus(varToRemove) == SCIP_VARSTATUS_AGGREGATED);
        assert(!infeasible && aggregated);
        SCIPfixVar(scip, edgeToVar({nodeToKeep, nodeToRemove}), 0.0, &infeasible, &fixed);
        assert(!infeasible && fixed);
    }

    void aggregateToDiffPartition(SCIP *scip, NetworKit::node nodeToRemove, NetworKit::node nodeToKeep,
                                  SCIP_Var *varToRemove, SCIP_Var *varToKeep) const {
        SCIP_Bool infeasible, redundant, aggregated, fixed;
        SCIPaggregateVars(scip, varToRemove, varToKeep, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated);
        assert(SCIPvarGetStatus(varToRemove) == SCIP_VARSTATUS_NEGATED);
        assert(!infeasible && aggregated);
        SCIPfixVar(scip, edgeToVar({nodeToKeep, nodeToRemove}), 1.0, &infeasible, &fixed);
        assert(!infeasible && fixed);
    }
};

/**
 * This function inits a SCIP
 * @param scip Output parameter
 * @param g The graph
 * @param root The sun
 */
void buildGeneralModel(Scip *s, const MaxCut &mc, bool allBinary);

ProbDataMc *getProbDataMc(SCIP *scip);

} // namespace sms

#endif // SMS_MC_EDGES_HPP
