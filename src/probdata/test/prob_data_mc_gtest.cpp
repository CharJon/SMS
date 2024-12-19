#include <gtest/gtest.h>
#include <scip/scip.h>

#include "sms/instance/maxcut.hpp"
#include "sms/probdata/prob_data_mc.hpp"

TEST(ProbDataMc, solveSquare) {
    auto g = NetworKit::Graph(2, true, false);
    g.addEdge(0, 1, -1);
    g.indexEdges();

    auto mcInstance = sms::MaxCut(g);

    // create scip with default plugins
    Scip *scip = nullptr;
    SCIP_CALL_EXC(SCIPcreate(&scip))
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip))
    ASSERT_NE(scip, nullptr);

    sms::buildGeneralModel(scip, mcInstance, true);

    auto solveReturnCode = SCIPsolve(scip);
    ASSERT_EQ(solveReturnCode, SCIP_OKAY);

    auto returnCode = SCIPfree(&scip);
    ASSERT_EQ(returnCode, SCIP_OKAY);
}
