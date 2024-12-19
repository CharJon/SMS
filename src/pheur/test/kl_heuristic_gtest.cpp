#include <filesystem>

#include <gtest/gtest.h>

#include "sms/io/io.hpp"
#include "sms/pheur/kl_heuristic.hpp"

TEST(KLHeuristic, ConstructorTest) {

    NetworKit::Graph g(4, true, false);

    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 1);

    sms::KLHeuristic kl(&g);
}

TEST(KLHeuristic, SimpleTest) {

    NetworKit::Graph g(4, true, false);

    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 2.0);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 1);

    sms::KLHeuristic kl(&g);

    std::vector<bool> allOnOneSide(g.upperNodeIdBound(), false);
    auto part = sms::BiPartition(allOnOneSide);
    auto dValues = kl.improvementValuesFromScratch(part);
    EXPECT_EQ(dValues[2], 3.0);

    kl.phase1Optimization(allOnOneSide, {2, 3});
    kl.phase2Optimization();
    auto sol = kl.getPrimalSolution();
    EXPECT_FALSE(sol[2]);
    EXPECT_FALSE(sol[3]);
    ASSERT_EQ(sol.size(), g.upperNodeIdBound());
}

TEST(KLHeuristic, SimpleRepartitionTest) {

    NetworKit::Graph g(4, true, false);

    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 1);

    sms::KLHeuristic kl(&g);

    kl.phase2Optimization();
    auto sol = kl.getPrimalSolution();
    ASSERT_EQ(sol.size(), g.upperNodeIdBound());
}

TEST(KLHeuristic, FixedVertices) {

    NetworKit::Graph g(6, true, false);

    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(2, 4, -1);
    g.addEdge(3, 4, 5);
    g.addEdge(0, 5, 1);
    g.addEdge(1, 5, 1);
    g.addEdge(4, 5, 1);

    sms::KLHeuristic kl(&g);

    std::vector<bool> init(g.upperNodeIdBound(), false);
    init[0] = true;
    kl.phase1Optimization(init, {0, 1});

    auto solution0 = kl.getPrimalSolution();

    EXPECT_TRUE(solution0[0]);
    EXPECT_TRUE(!solution0[1]);
    kl.phase2Optimization();
    auto solution1 = kl.getPrimalSolution();
    EXPECT_TRUE(solution0[0]);
    EXPECT_TRUE(!solution0[1]);

    ASSERT_EQ(solution0.size(), g.upperNodeIdBound());
    ASSERT_EQ(solution1.size(), g.upperNodeIdBound());
}