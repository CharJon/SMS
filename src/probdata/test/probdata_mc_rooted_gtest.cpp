#include <gtest/gtest.h>
#include <scip/scip.h>

#include "sms/graph/graphs.hpp"
#include "sms/io/io.hpp"
#include "sms/probdata/mc_rooted.hpp"

TEST(ProbDataRootedMc, solvePm1s_50_1) {
    NetworKit::Graph g = mcFileToGraph("test/data/pm1s_50.1.mc");
    auto r = sms::compactGraph(g);
    g = std::move(r.compactGraph);
    g.indexEdges();

    Scip *scip = nullptr;
    buildRootMcModel(&scip, &g, 1, false);
    ASSERT_NE(scip, nullptr);
    SCIPsolve(scip);
    SCIPfree(&scip);
}

TEST(ProbDataRootedMc, solveSquare) {
    sleep(1);
    NetworKit::Graph g = mcFileToGraph("test/data/square.mc");
    g.indexEdges();

    Scip *scip = nullptr;
    buildRootMcModel(&scip, &g, 0, false);
    ASSERT_NE(scip, nullptr);
    SCIPsolve(scip);
    SCIPfree(&scip);
}

TEST(ProbDataRootedMc, solvePm1s_50_1_plugins) {
    NetworKit::Graph g = mcFileToGraph("test/data/pm1s_50.1.mc");
    auto r = sms::compactGraph(g);
    g = std::move(r.compactGraph);
    g.indexEdges();

    Scip *scip = nullptr;
    buildRootMcModel(&scip, &g, 1, false);
    ASSERT_NE(scip, nullptr);

    SCIPsolve(scip);

    SCIPfree(&scip);
}
