#include <gtest/gtest.h>

#include "networkit/io/EdgeListReader.hpp"

#include "sms/presol/presolver_mc.hpp"

TEST(Presolver, Star) {
    // tree has only bridges as biconnected components
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, -1);
    g.addEdge(0, 3, -1);

    sms::PresolverMC presolver(g);
    presolver.run();
    EXPECT_EQ(presolver.subgraphs().size(), 3);
    for (const auto &subgraph : presolver.subgraphs()) {
        EXPECT_EQ(subgraph.numberOfNodes(), 2);
    }

    std::vector<sms::Partition> emptySol(2, sms::Partition::kUNASSIGNED);
    presolver.addSolution(0, emptySol);
    presolver.addSolution(1, emptySol);
    presolver.addSolution(2, emptySol);

    auto fullSol = presolver.recoverSolution();
    for (auto s : fullSol) {
        EXPECT_EQ(s, false);
    }
}

TEST(Presolver, Disconnected) {
    // tree has only bridges as biconnected components
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(2, 3, -1);

    sms::PresolverMC presolver(g);
    presolver.run();
    EXPECT_EQ(presolver.subgraphs().size(), 2);
}

TEST(Presolver, Bigger) {
    // tree has only bridges as biconnected components
    NetworKit::Graph g(10, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, -1);
    g.addEdge(1, 2, -1);
    g.addEdge(1, 3, -1);
    g.addEdge(1, 4, -1);
    g.addEdge(2, 5, -1);
    g.addEdge(2, 6, -1);
    g.addEdge(2, 7, -1);
    g.addEdge(3, 4, -1);
    g.addEdge(5, 6, -1);
    g.addEdge(7, 8, -1);
    g.addEdge(7, 9, -1);
    g.addEdge(8, 9, -1);
    ASSERT_EQ(g.numberOfNodes(), 10);
    ASSERT_EQ(g.numberOfEdges(), 13);

    sms::PresolverMC presolver(g);
    presolver.run();
    EXPECT_EQ(presolver.subgraphs().size(), 5);
    for (unsigned int i = 0; i < presolver.subgraphs().size(); i++) {
        const auto root = presolver.roots()[i];
        const auto &subgraph = presolver.subgraphs()[i];
        // const auto &sub2orig = presolver.toOrig()[i];
        bool rootFound = false;
        for (auto u : subgraph.nodeRange()) {
            rootFound |= (u == root);
        }
        EXPECT_TRUE(rootFound);
    }

    // when traversing subgraphs in order, roots always have a set value
    std::vector<sms::Partition> sols(g.upperNodeIdBound(), sms::Partition::kUNASSIGNED);
    for (auto u : presolver.subgraphs()[0].nodeRange()) {
        sols[presolver.toOrig()[0][u]] = sms::Partition::kZERO;
    }

    for (unsigned int i = 1; i < presolver.subgraphs().size(); i++) {
        const auto &subgraph = presolver.subgraphs()[i];
        const auto root = presolver.roots()[i];
        const auto &sub2orig = presolver.toOrig()[i];
        EXPECT_NE(sols[root], sms::Partition::kUNASSIGNED);
        for (auto u : subgraph.nodeRange()) {
            sols[sub2orig[u]] = sms::Partition::kONE;
        }
    }
}
