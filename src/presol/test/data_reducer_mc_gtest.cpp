#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "networkit/io/EdgeListReader.hpp"
#include "tlx/container/d_ary_addressable_int_heap.hpp"

#include "sms/io/io.hpp"
#include "sms/presol/data_reducer_mc.hpp"
#include "sms/solver/enumeration_solver.hpp"

NetworKit::edgeweight bestSolutionFromNothing(const NetworKit::Graph &g, const sms::DataReducerMC &dr) {
    auto sol = std::vector<sms::Partition>(g.upperNodeIdBound(), sms::Partition::kUNASSIGNED);
    auto origSol = dr.recoverSolution(sol);
    auto origSolValue = sms::solutionValue(g, origSol);
    return origSolValue;
}

TEST(DataReduction, RuleFlags) {
    std::cout << sms::to_underlying_type(sms::RuleFlag::kAll) << std::endl;
    EXPECT_TRUE(sms::RuleFlag::kAll >= sms::RuleFlag::kDegreeZero);
}

TEST(DataReduction, OneEdge) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, -1);

    sms::DataReducerMC dr(g);
    dr.deactivateAllRules(0);
    dr.deactivateAllRules(1);
    dr.deactivateAllRules(2);
    dr.activateRule(0, sms::RuleFlag::kDominatingEdge);
    ASSERT_EQ(dr.getTop(), 0);
    dr.activateRule(1, sms::RuleFlag::kDegreeZero);
    ASSERT_EQ(dr.getTop(), 1);
    auto next = dr.extractTop();
    ASSERT_EQ(next, 1);
    ASSERT_EQ(dr.getTop(), 0);
}

TEST(DataReduction, Line1) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 3, -1);
    g.addEdge(3, 4, 15);

    sms::DataReducerMC dr(g);
    dr.initRun();
    auto numAggregated = dr.nextNodeAggregation();
    ASSERT_GT(numAggregated, 0);
}

TEST(DataReduction, Line2) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 3, -1);
    g.addEdge(3, 4, 15);

    sms::DataReducerMC dr(g);
    dr.run();
    ASSERT_EQ(dr.getGraph().numberOfNodes(), 0);
    EXPECT_EQ(dr.getOffset(), 18);

    auto sol = std::vector<sms::Partition>(g.upperNodeIdBound(), sms::Partition::kUNASSIGNED);
    auto origSol = dr.recoverSolution(sol);
    bool flip = origSol[0] != sms::Partition::kZERO;
    EXPECT_TRUE(origSol[1] == (sms::Partition::kZERO ^ flip));
    EXPECT_TRUE(origSol[2] == (sms::Partition::kONE ^ flip));
    EXPECT_TRUE(origSol[3] == (sms::Partition::kONE ^ flip));
    EXPECT_TRUE(origSol[4] == (sms::Partition::kZERO ^ flip));
}

TEST(DataReduction, Cycle1) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 3, -1);
    g.addEdge(3, 4, 15);
    g.addEdge(4, 1, 7);

    sms::DataReducerMC dr(g);
    dr.initRun();
    auto numAggregated = dr.nextNodeAggregation();
    ASSERT_GT(numAggregated, 0);
}

TEST(DataReduction, Cycle2) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 3, -1);
    g.addEdge(3, 4, 15);
    g.addEdge(4, 1, 7);

    sms::DataReducerMC dr(g);
    dr.run();
    ASSERT_EQ(dr.getGraph().numberOfNodes(), 0);
    EXPECT_EQ(dr.getOffset(), 24);

    auto sol = std::vector<sms::Partition>(g.upperNodeIdBound(), sms::Partition::kUNASSIGNED);
    auto origSol = dr.recoverSolution(sol);
    bool flip = origSol[0] != sms::Partition::kZERO;
    EXPECT_TRUE(origSol[1] == (sms::Partition::kZERO ^ flip));
    EXPECT_TRUE(origSol[2] == (sms::Partition::kONE ^ flip));
    EXPECT_TRUE(origSol[3] == (sms::Partition::kZERO ^ flip));
    EXPECT_TRUE(origSol[4] == (sms::Partition::kONE ^ flip));
}

TEST(DataReduction, TwoMergedTriangles) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 0, -1);
    g.addEdge(1, 3, 7);
    g.addEdge(3, 2, -7);

    sms::DataReducerMC dr(g);
    dr.initRun();
    auto numAggregated = dr.nextNodeAggregation();
    ASSERT_GT(numAggregated, 0);
    dr.run();
    EXPECT_EQ(dr.getOffset(), 9);
}

TEST(DataReduction, FourNodes) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, -1);
    g.addEdge(0, 3, -1);
    g.addEdge(1, 2, -1);

    sms::DataReducerMC dr(g);
    dr.initRun();
    auto numAggregated = dr.nextNodeAggregation();
    ASSERT_GT(numAggregated, 0);
    numAggregated = dr.nextNodeAggregation();
    ASSERT_GT(numAggregated, 0);
}

TEST(DataReduction, DegreeThreePlus1) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);

    sms::DataReducerMC dr(g, 3);
    auto numAggregated = dr.degreeThreeApplication(0);
    ASSERT_GT(numAggregated, 0);
    // ASSERT_FALSE(dr.getHistory()[0]);
    EXPECT_EQ(dr.getOffset(), 3);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 2));
    EXPECT_EQ(dr.getGraph().weight(1, 2), -0.5);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 3));
    EXPECT_EQ(dr.getGraph().weight(1, 3), -0.5);
    EXPECT_TRUE(dr.getGraph().hasEdge(2, 3));
    EXPECT_EQ(dr.getGraph().weight(2, 3), -0.5);

    auto dr2 = sms::DataReducerMC(g);
    dr2.run();
    auto solOrig = bestSolutionFromNothing(g, dr2);
    EXPECT_EQ(solOrig, dr2.getOffset());
}

TEST(DataReduction, DegreeThreeMinus1) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, -1);
    g.addEdge(0, 3, -1);

    sms::DataReducerMC dr(g, 3);
    auto numAggregated = dr.degreeThreeApplication(0);
    ASSERT_GT(numAggregated, 0);
    // ASSERT_FALSE(dr.getHistory()[0].simple);
    EXPECT_EQ(dr.getOffset(), 0);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 2));
    EXPECT_EQ(dr.getGraph().weight(1, 2), -0.5);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 3));
    EXPECT_EQ(dr.getGraph().weight(1, 3), -0.5);
    EXPECT_TRUE(dr.getGraph().hasEdge(2, 3));
    EXPECT_EQ(dr.getGraph().weight(2, 3), -0.5);
}

TEST(DataReduction, DegreeThreeMixed1) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 5);
    g.addEdge(0, 2, 6);
    g.addEdge(0, 3, 7);

    sms::DataReducerMC dr(g, 3);
    auto numAggregated = dr.degreeThreeApplication(0);
    ASSERT_GT(numAggregated, 0);
    // ASSERT_FALSE(dr.getHistory()[0].simple);
    EXPECT_EQ(dr.getOffset(), 18);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 2));
    EXPECT_EQ(dr.getGraph().weight(1, 2), -2.0);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 3));
    EXPECT_EQ(dr.getGraph().weight(1, 3), -3.0);
    EXPECT_TRUE(dr.getGraph().hasEdge(2, 3));
    EXPECT_EQ(dr.getGraph().weight(2, 3), -4.0);
}

TEST(DataReduction, DegreeThreeMixed2) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 5);
    g.addEdge(0, 2, -6);
    g.addEdge(0, 3, 7);

    sms::DataReducerMC dr(g, 3);
    auto numAggregated = dr.degreeThreeApplication(0);
    ASSERT_GT(numAggregated, 0);
    // ASSERT_FALSE(dr.getHistory()[0].simple);
    ASSERT_EQ(dr.getHistory()[0].toAggregate, 0);
    EXPECT_EQ(dr.getOffset(), 6);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 2));
    EXPECT_EQ(dr.getGraph().weight(1, 2), 2.0);
    EXPECT_TRUE(dr.getGraph().hasEdge(1, 3));
    EXPECT_EQ(dr.getGraph().weight(1, 3), -3.0);
    EXPECT_TRUE(dr.getGraph().hasEdge(2, 3));
    EXPECT_EQ(dr.getGraph().weight(2, 3), 4.0);
}

TEST(DataReduction, UnitWeigthCycle) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 0, 1);

    sms::DataReducerMC dr(g);
    dr.run();
    ASSERT_EQ(dr.getGraph().numberOfNodes(), 0);
    EXPECT_EQ(dr.getOffset(), 4);
    auto solOrig = bestSolutionFromNothing(g, dr);
    EXPECT_EQ(solOrig, dr.getOffset());
}

TEST(DataReduction, ThreeRegular) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 5);
    g.addEdge(0, 2, -6);
    g.addEdge(0, 3, 7);
    g.addEdge(1, 2, -7);
    g.addEdge(1, 3, 8);
    g.addEdge(2, 3, -9);

    sms::DataReducerMC dr(g, 3);
    dr.run();
    ASSERT_EQ(dr.getGraph().numberOfNodes(), 0);
    EXPECT_EQ(dr.getOffset(), 6);
    auto solOrig = bestSolutionFromNothing(g, dr);
    EXPECT_EQ(solOrig, dr.getOffset());
}

TEST(DataReduction, NodeMapping) {
    NetworKit::Graph g(0, true, false);

    sms::DataReducerMC dr(g);

    std::array<NetworKit::node, 3> neighbors = {1, 2, 3};

    std::array<uint8_t, 4> column = {0, 0, 0, 0};
    auto res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 1);
    EXPECT_TRUE(res.reverse({true}));
    EXPECT_FALSE(res.reverse({false}));

    column = {1, 1, 1, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 1);
    EXPECT_TRUE(res.reverse({false}));
    EXPECT_FALSE(res.reverse({true}));

    column = {0, 0, 1, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 2);
    EXPECT_TRUE(res.reverse({true}));
    EXPECT_FALSE(res.reverse({false}));

    column = {1, 1, 0, 0};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 2);
    EXPECT_TRUE(res.reverse({false}));
    EXPECT_FALSE(res.reverse({true}));

    column = {0, 1, 0, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 3);
    EXPECT_TRUE(res.reverse({true}));
    EXPECT_FALSE(res.reverse({false}));

    column = {1, 0, 1, 0};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    EXPECT_EQ(res.others[0], 3);
    EXPECT_TRUE(res.reverse({false}));
    EXPECT_FALSE(res.reverse({true}));

    column = {0, 1, 0, 0};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_FALSE(res.reverse({false, false, false}));
    EXPECT_TRUE(res.reverse({false, false, true}));
    EXPECT_FALSE(res.reverse({false, true, false}));
    EXPECT_FALSE(res.reverse({false, true, true}));

    EXPECT_TRUE(res.reverse({true, true, true}));
    EXPECT_FALSE(res.reverse({true, true, false}));
    EXPECT_TRUE(res.reverse({true, false, true}));
    EXPECT_TRUE(res.reverse({true, false, false}));

    column = {1, 0, 1, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_TRUE(res.reverse({false, false, false}));
    EXPECT_FALSE(res.reverse({false, false, true}));
    EXPECT_TRUE(res.reverse({false, true, false}));
    EXPECT_TRUE(res.reverse({false, true, true}));

    EXPECT_FALSE(res.reverse({true, true, true}));
    EXPECT_TRUE(res.reverse({true, true, false}));
    EXPECT_FALSE(res.reverse({true, false, true}));
    EXPECT_FALSE(res.reverse({true, false, false}));

    column = {0, 0, 1, 0};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_FALSE(res.reverse({false, false, false}));
    EXPECT_FALSE(res.reverse({false, false, true}));
    EXPECT_TRUE(res.reverse({false, true, false}));
    EXPECT_FALSE(res.reverse({false, true, true}));

    EXPECT_TRUE(res.reverse({true, true, true}));
    EXPECT_TRUE(res.reverse({true, true, false}));
    EXPECT_FALSE(res.reverse({true, false, true}));
    EXPECT_TRUE(res.reverse({true, false, false}));

    column = {1, 1, 0, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_TRUE(res.reverse({false, false, false}));
    EXPECT_TRUE(res.reverse({false, false, true}));
    EXPECT_FALSE(res.reverse({false, true, false}));
    EXPECT_TRUE(res.reverse({false, true, true}));

    EXPECT_FALSE(res.reverse({true, true, true}));
    EXPECT_FALSE(res.reverse({true, true, false}));
    EXPECT_TRUE(res.reverse({true, false, true}));
    EXPECT_FALSE(res.reverse({true, false, false}));

    column = {0, 0, 0, 1};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_FALSE(res.reverse({false, false, false}));
    EXPECT_FALSE(res.reverse({false, false, true}));
    EXPECT_FALSE(res.reverse({false, true, false}));
    EXPECT_TRUE(res.reverse({false, true, true}));

    EXPECT_TRUE(res.reverse({true, true, true}));
    EXPECT_TRUE(res.reverse({true, true, false}));
    EXPECT_TRUE(res.reverse({true, false, true}));
    EXPECT_FALSE(res.reverse({true, false, false}));

    column = {1, 1, 1, 0};
    res = sms::columnToNodeAggregation(0, neighbors, column);
    // ASSERT_FALSE(res.simple);

    EXPECT_TRUE(res.reverse({false, false, false}));
    EXPECT_TRUE(res.reverse({false, false, true}));
    EXPECT_TRUE(res.reverse({false, true, false}));
    EXPECT_FALSE(res.reverse({false, true, true}));

    EXPECT_FALSE(res.reverse({true, true, true}));
    EXPECT_FALSE(res.reverse({true, true, false}));
    EXPECT_FALSE(res.reverse({true, false, true}));
    EXPECT_TRUE(res.reverse({true, false, false}));
}

TEST(DataReduction, egoFacebook) {
    NetworKit::EdgeListReader el(' ', 1, "#", true);
    auto g = el.read("test/data/ego-facebook.wel");
    sms::DataReducerMC dr(g);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 12); // trivial value for removing nodes with degree <= 2
    if (dr.getGraph().numberOfNodes() == 0) {
        EXPECT_EQ(dr.getOffset(), 2975);
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_EQ(origSolValue, 2975);
    }
}

TEST(DataReduction, imgseg_271031) {
    NetworKit::Graph g = mcFileToGraph("test/data/imgseg_271031.txt");
    g.restoreNode(0); // fastest way to make g compact

    sms::DataReducerMC dr(g, 3);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 0); // all degrees plus dominating edge
    if (dr.getGraph().numberOfNodes() == 0) {
        double optimalSolution = 3016231197928 / 1e6;
        double tolerance = 1e-5;
        EXPECT_NEAR(dr.getOffset(), optimalSolution, tolerance);
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_NEAR(origSolValue, optimalSolution, tolerance);
    }
}

TEST(DataReduction, imgseg_35058) {
    NetworKit::Graph g = mcFileToGraph("test/data/imgseg_35058.txt");
    g.restoreNode(0); // fastest way to make g compact

    sms::DataReducerMC dr(g, 3);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 17); // all degrees plus dominating edge
    if (dr.getGraph().numberOfNodes() == 0) {
        double optimalSolution = 1331964975484 / 1e6;
        double tolerance = 1e-5;
        EXPECT_NEAR(dr.getOffset(), optimalSolution, tolerance);
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_NEAR(origSolValue, optimalSolution, tolerance);
    }
}

TEST(DataReduction, imgseg_106025) {
    NetworKit::Graph g = mcFileToGraph("test/data/imgseg_106025.txt");
    g.restoreNode(0); // fastest way to make g compact

    sms::DataReducerMC dr(g, 3);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 21); // all degrees plus dominating edge
    if (dr.getGraph().numberOfNodes() == 0) {
        double optimalSolution = 2927960456661 / 1e6;
        double tolerance = 1e-4;
        EXPECT_NEAR(dr.getOffset(), optimalSolution, tolerance);
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_NEAR(origSolValue, optimalSolution, tolerance);
    }
}

TEST(DataReduction, imgseg_374020) {
    NetworKit::Graph g = mcFileToGraph("test/data/imgseg_374020.mc");
    g.restoreNode(0); // fastest way to make g compact

    sms::DataReducerMC dr(g, 3);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 109); // without degree three
    if (dr.getGraph().numberOfNodes() == 0) {
        ASSERT_EQ(dr.getHistory().size(), g.numberOfNodes());
        double optimalSolution = 8074339398623 / 1e6;
        double tolerance = 1e-4;
        EXPECT_NEAR(dr.getOffset(), optimalSolution, tolerance);
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_NEAR(origSolValue, dr.getOffset(), tolerance);
    }
}

TEST(DataReduction, imgseg_105019) {
    NetworKit::Graph g = mcFileToGraph("test/data/imgseg_105019.txt");
    g.restoreNode(0); // fastest way to make g compact

    sms::DataReducerMC dr(g, 3);
    dr.run();
    EXPECT_LE(dr.getGraph().numberOfNodes(), 0); // all degrees plus dominating edge
    if (dr.getGraph().numberOfNodes() == 0) {
        double optimalSolution = 7381459364666 / 1e6;
        double tolerance = 1e-4;
        EXPECT_NEAR(dr.getOffset(), optimalSolution, tolerance);
        // nothing left, so can recover from nothing
        auto origSolValue = bestSolutionFromNothing(g, dr);
        EXPECT_NEAR(origSolValue, optimalSolution, tolerance);
    }
}

TEST(DataReduction, reverse) {
    auto g = NetworKit::Graph(3, true, false);

    auto dr = sms::DataReducerMC(g);

    std::vector<sms::NodeAggregation> aggr;
    aggr.emplace_back(0, 1, false);
    aggr.emplace_back(1, 2, true);
    aggr.emplace_back(2, 3, false);

    std::vector<sms::Partition> transformedSol = {sms::Partition::kUNASSIGNED, sms::Partition::kUNASSIGNED,
                                                  sms::Partition::kUNASSIGNED, sms::Partition::kZERO};

    auto res = sms::recoverSolution(aggr, transformedSol);

    EXPECT_THAT(res, testing::ElementsAre(sms::Partition::kONE, sms::Partition::kONE, sms::Partition::kZERO,
                                          sms::Partition::kZERO));
}

TEST(DataReduction, GomoryHu) {
    auto g = NetworKit::Graph(4, true, false, true);
    g.addEdge(0, 1, -1);
    g.addEdge(1, 3, 2);
    g.addEdge(0, 2, 1);
    g.addEdge(2, 3, 2);

    auto dr = sms::DataReducerMC(g);

    dr.gomoryHuApplication();

    auto compactMapId = NetworKit::GraphTools::getContinuousNodeIds(dr.getGraph());
    auto compactG = NetworKit::GraphTools::getCompactedGraph(dr.getGraph(), compactMapId);
    sms::EnumerationSolver solver(compactG);

    solver.run();
    auto maxCut = solver.bestSolutionValue();

    EXPECT_EQ(4, maxCut + dr.getOffset());
}

TEST(DataReduction, SeparatedClique1) {
    auto g = NetworKit::Graph(7, true, false, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);

    for (unsigned int i = 4; i < 7; i++) {
        for (unsigned int j = 0; j < i; j++) {
            g.addEdge(i, j, 1);
        }
    }

    auto dr = sms::DataReducerMC(g);
    dr.initRun();
    dr.removeSeparatedCliques();
    EXPECT_EQ(dr.getGraph().numberOfNodes(), 4);
    EXPECT_EQ(dr.getHistory().size(), 3);
    for (const auto &contraction : dr.getHistory()) {
        // EXPECT_TRUE(contraction.simple);
    }
    EXPECT_FALSE(dr.getGraph().hasEdge(0, 1));
    ASSERT_TRUE(dr.getGraph().hasEdge(1, 3));
    ASSERT_EQ(dr.getGraph().weight(1, 3), -1);
}
