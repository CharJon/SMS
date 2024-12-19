#ifndef SMS_DISP_HPP
#define SMS_DISP_HPP

#include <iostream>

#include "objscip/objscip.h"

namespace sms {

class ExtendedDisp : public scip::ObjDisp {

public:
    /** constructor */
    explicit ExtendedDisp(Scip *scip);

    /** destructor */
    ~ExtendedDisp() override = default;

    /** destructor of display column to free user data (called when SCIP is exiting)
     *
     *  @see SCIP_DECL_DISPFREE(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_free(SCIP * /*scip*/, SCIP_DISP * /*disp*/) override { return SCIP_OKAY; }

    /** initialization method of display column (called after problem was transformed)
     *
     *  @see SCIP_DECL_DISPINIT(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_init(SCIP * /*scip*/, SCIP_DISP * /*disp*/) override { return SCIP_OKAY; }

    /** deinitialization method of display column (called before transformed problem is freed)
     *
     *  @see SCIP_DECL_DISPEXIT(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_exit(SCIP * /*scip*/, SCIP_DISP * /*disp*/) override { return SCIP_OKAY; }

    /** solving process initialization method of display column (called when branch and bound process is about to begin)
     *
     *  @see SCIP_DECL_DISPINITSOL(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_initsol(SCIP * /*scip*/, SCIP_DISP * /*disp*/) override { return SCIP_OKAY; }

    /** solving process deinitialization method of display column (called before branch and bound process data is freed)
     *
     *  @see SCIP_DECL_DISPEXITSOL(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_exitsol(SCIP * /*scip*/, SCIP_DISP * /*disp*/) override { return SCIP_OKAY; }

    /** output method of display column to output file stream 'file'
     *
     *  @see SCIP_DECL_DISPOUTPUT(x) in @ref type_disp.h
     */
    SCIP_RETCODE scip_output(SCIP *scip, SCIP_DISP * /*disp*/, FILE *file) override;
};

} // namespace sms

#endif // SMS_DISP_HPP
