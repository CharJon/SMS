#include <gtest/gtest.h>

#include "sms/graph/gomory_hu_tree.hpp"

TEST(GomoryHuTree, BasicTree) {
    auto g = NetworKit::Graph(6, true, false, true);
    g.addEdge(0, 1, 7.3);
    g.addEdge(0, 2, 7.0);
    g.addEdge(1, 3, 10.0);
    g.addEdge(3, 4, 1.0);
    g.addEdge(3, 5, 5.0);
    g.indexEdges();
    auto ght = sms::GomoryHuTree(g);
    ght.run();
    NetworKit::Graph gomoryHuTree = ght.getGomoryHuTree();
    // Check edges
    ASSERT_EQ(ght.minCutEdge(0, 1).weight, 7.3);
    ASSERT_EQ(ght.minCutEdge(0, 2).weight, 7.0);
    ASSERT_EQ(ght.minCutEdge(1, 3).weight, 10.0);
    ASSERT_EQ(ght.minCutEdge(3, 4).weight, 1.0);
    ASSERT_EQ(ght.minCutEdge(3, 5).weight, 5.0);
    // check not adjacent edge pairs
    ASSERT_EQ(ght.minCutEdge(2, 4).weight, 1.0);
    ASSERT_EQ(ght.minCutEdge(0, 5).weight, 5.0);

    // check that gomoryHuTree has the right amount of numbers and nodes
    ASSERT_EQ(gomoryHuTree.numberOfNodes(), 6);
    ASSERT_EQ(gomoryHuTree.numberOfEdges(), 5);
}

TEST(GomoryHuTree, BasicCircle) {
    auto g = NetworKit::Graph(4, true, false, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 3, 2);
    g.addEdge(2, 3, 3);
    g.addEdge(0, 2, 4);
    auto ght = sms::GomoryHuTree(g);
    ght.run();
    NetworKit::Graph gomoryHuTree = ght.getGomoryHuTree();
    ASSERT_EQ(ght.minCutEdge(0, 1).weight, 3);
    ASSERT_EQ(ght.minCutEdge(0, 2).weight, 5);
    ASSERT_EQ(ght.minCutEdge(2, 3).weight, 4);
}

TEST(GomoryHuTree, FullyConnectedGraph) {
    auto g = NetworKit::Graph(4, true, false, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 3, 2);
    g.addEdge(2, 3, 3);
    g.addEdge(0, 2, 4);
    g.addEdge(0, 3, 5);
    g.addEdge(1, 2, 5);
    auto ght = sms::GomoryHuTree(g);
    ght.run();
    NetworKit::Graph gomoryHuTree = ght.getGomoryHuTree();
    ASSERT_EQ(ght.minCutEdge(0, 1).weight, 8);
    ASSERT_EQ(ght.minCutEdge(0, 2).weight, 10);
    ASSERT_EQ(ght.minCutEdge(2, 3).weight, 10);
}

TEST(GomoryHuTree, TwoComponents) {
    auto g = NetworKit::Graph(9, true, false, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 3, 2);
    g.addEdge(2, 3, 3);
    g.addEdge(0, 2, 4);
    g.addEdge(0, 3, 5);
    g.addEdge(1, 2, 5);

    g.addEdge(7, 5, 1);
    g.addEdge(7, 6, 2);

    auto ght = sms::GomoryHuTree(g);
    ght.run();
    NetworKit::Graph gomoryHuTree = ght.getGomoryHuTree();
    ASSERT_EQ(ght.minCutEdge(0, 1).weight, 8);
    ASSERT_EQ(ght.minCutEdge(0, 2).weight, 10);
    ASSERT_EQ(ght.minCutEdge(2, 3).weight, 10);

    ASSERT_EQ(ght.minCutEdge(5, 6).weight, 1);
    ASSERT_EQ(ght.minCutEdge(2, 5).weight, 0);
    ASSERT_EQ(ght.minCutEdge(2, 4).weight, 0);
}