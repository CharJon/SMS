#ifndef SMS_PHEUR_KL_HPP
#define SMS_PHEUR_KL_HPP

#include "objscip/objscip.h"
#include "scip/scip.h"

#include "networkit/graph/Graph.hpp"

#include "sms/pheur/kl_heuristic.hpp"

namespace sms {

class PHeurKl : public scip::ObjHeur {

private:
    NetworKit::Graph const *const mcGraph_;
    KLHeuristic heuristic_;
    int lastSolutionIndex_ = 0; // keep track of which solutions from the solution pool have been checked already

public:
    /** default constructor */
    explicit PHeurKl(SCIP *scip, NetworKit::Graph const *g);

    /** destructor */
    ~PHeurKl() override = default;

    SCIP_RETCODE improveSolution(SCIP *scip, SCIP_HEUR *heur, SCIP_SOL *sol, SCIP_RESULT *result);

    /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
    SCIP_DECL_HEURFREE(scip_free) override;

    /** initialization method of primal heuristic (called after problem was transformed) */
    SCIP_DECL_HEURINIT(scip_init) override;

    /** deinitialization method of primal heuristic (called before transformed problem is freed) */
    SCIP_DECL_HEUREXIT(scip_exit) override;

    /** solving process initialization method of primal heuristic (called when branch and bound process is about to
     * begin)
     *
     *  This method is called when the presolving was finished and the branch and bound process is about to begin.
     *  The primal heuristic may use this call to initialize its branch and bound specific data.
     *
     */
    SCIP_DECL_HEURINITSOL(scip_initsol) override;

    /** solving process deinitialization method of primal heuristic (called before branch and bound process data is
     * freed)
     *
     *  This method is called before the branch and bound process is freed.
     *  The primal heuristic should use this call to clean up its branch and bound data.
     */
    SCIP_DECL_HEUREXITSOL(scip_exitsol) override;

    /** execution method of primal heuristic
     *
     *  Searches for feasible primal solutions. The method is called in the node processing loop.
     *
     *  possible return values for *result:
     *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
     *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
     *  - SCIP_DIDNOTRUN  : the heuristic was skipped
     *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
     *                      its frequency
     */
    SCIP_DECL_HEUREXEC(scip_exec) override;

    /** clone method which will be used to copy a objective plugin */
    SCIP_DECL_HEURCLONE(ObjCloneable *clone) override; /*lint !e665*/

    /** returns whether the objective plugin is copyable */
    SCIP_DECL_HEURISCLONEABLE(iscloneable) override { return false; }

}; /*lint !e1712*/

} // namespace sms

#endif // SMS_PHEUR_KL_HPP
