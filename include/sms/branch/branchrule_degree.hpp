#ifndef SMS_BRANCHRULEDEGREE_HPP
#define SMS_BRANCHRULEDEGREE_HPP

#include "objscip/objscip.h"

class BranchruleDegree : public scip::ObjBranchrule {

public:
    explicit BranchruleDegree(SCIP *scip);

    SCIP_RETCODE
    scip_execlp(SCIP *scip, SCIP_BRANCHRULE *branchrule, unsigned int allowaddcons, SCIP_RESULT *result) override;

private:
    SCIP_Bool dynamic_;
    SCIP_Bool tieBreak_;

    int degree(SCIP_VAR *variable);
};

#endif // SMS_BRANCHRULEDEGREE_HPP
