#include "sms/disp/disp.hpp"

#define DISP_NAME "ocw"
#define DISP_DESC "odd-cycle"
#define DISP_HEADER "cuts"
#define DISP_WIDTH 21       /**< the width of the display column */
#define DISP_PRIORITY 10000 /**< the priority of the display column */
#define DISP_POSITION 30300 /**< the relative position of the display column */
#define DISP_STRIPLINE                                                                                                 \
    TRUE /**< the default for whether the display column should be separated with a line from its right neighbor */

namespace sms {

ExtendedDisp::ExtendedDisp(Scip *scip)
    : scip::ObjDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) {
}

SCIP_RETCODE ExtendedDisp::scip_output(SCIP *scip, SCIP_DISP *, FILE *file) {
    // assert(disp != NULL);
    assert(scip != nullptr);

    // MightDo: Use SCIPinfoMessage() instead of SCIPmessageFPrintInfo()

    auto triConsHndlr = SCIPfindConshdlr(scip, "triangleparity");
    if (triConsHndlr != nullptr) {
        auto totalNumCutsFound = SCIPconshdlrGetNCutsFound(triConsHndlr);
        SCIPdispLongint(SCIPgetMessagehdlr(scip), file, totalNumCutsFound, 4);
        SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " | ");
    }
    auto holesConsHndlr = SCIPfindConshdlr(scip, "holeparity");
    if (holesConsHndlr != nullptr) {
        auto totalNumCutsFound = SCIPconshdlrGetNCutsFound(holesConsHndlr);
        SCIPdispLongint(SCIPgetMessagehdlr(scip), file, totalNumCutsFound, 4);
        SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " | ");
    }
    auto cliqueConsHndlr = SCIPfindConshdlr(scip, "cliquecut");
    if (cliqueConsHndlr != nullptr) {
        auto totalNumCutsFound = SCIPconshdlrGetNCutsFound(cliqueConsHndlr);
        SCIPdispLongint(SCIPgetMessagehdlr(scip), file, totalNumCutsFound, 4);
        SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " | ");
    }

    return SCIP_OKAY;
}

} // namespace sms
