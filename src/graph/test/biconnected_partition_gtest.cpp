#include <gtest/gtest.h>

#include "sms/graph/biconnected_partition.hpp"
#include "sms/io/io.hpp"

TEST(BiconnectedPartition, SimpleGraph1) {
    auto g = NetworKit::Graph(5);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(2, 4);
    g.removeNode(3);

    sms::BiconnectedPartition biconn(g);
    biconn.run();

    EXPECT_EQ(biconn.components().size(), 3);
    EXPECT_EQ(biconn.numComponents(), 3);
    EXPECT_EQ(biconn.articulationPoints().size(), 2);
    EXPECT_EQ(biconn.numArticulationPoints(), 2);
}

TEST(BiconnectedPartition, SimpleGraph2) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(0, 4);
    g.addEdge(3, 4);
    g.addEdge(1, 2);

    sms::BiconnectedPartition biconn(g);
    biconn.run();

    ASSERT_EQ(biconn.numComponents(), 2);
    ASSERT_EQ(biconn.numArticulationPoints(), 1);

    ASSERT_EQ(biconn.tree().numberOfNodes(), 3);
    ASSERT_EQ(biconn.tree().numberOfEdges(), 2);

    for (unsigned int i = 0; i < biconn.numComponents(); i++) {
        EXPECT_FALSE(biconn.isArticulationPoint(i));
    }

    for (unsigned int i = biconn.numComponents(); i < biconn.numComponents() + biconn.numArticulationPoints(); i++) {
        EXPECT_TRUE(biconn.isArticulationPoint(i));
    }
}

TEST(BiconnectedPartition, GraphWithBridge) {
    NetworKit::Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);

    g.addEdge(0, 5);

    g.addEdge(5, 4);
    g.addEdge(3, 4);
    g.addEdge(3, 5);

    sms::BiconnectedPartition biconn(g);
    biconn.run();

    ASSERT_EQ(biconn.components().size(), 3);
    // ASSERT_EQ(biconn.bridges().size(), 1);

    // ASSERT_EQ(biconn.tree.numberOfNodes(), 5);
    // ASSERT_EQ(biconn.tree.numberOfEdges(), 4);
    // ASSERT_EQ(biconn.tree_map[0].size(), 3);

    // ASSERT_EQ(biconn.inverse_tree_map[0], 3);
    // ASSERT_EQ(biconn.inverse_tree_map[5], 4);
}

TEST(BiconnectedPartition, StarGraph) {
    NetworKit::Graph g = mcFileToGraph("test/data/star.mc");

    sms::BiconnectedPartition biconn(g);
    biconn.run();

    ASSERT_EQ(biconn.components().size(), 20);
    // ASSERT_EQ(biconn.bridges().size(), 20);

    // ASSERT_EQ(biconn.tree.numberOfNodes(), 21);
    // ASSERT_EQ(biconn.tree.numberOfEdges(), 20);
    // ASSERT_EQ(biconn.tree_map[0].size(), 2);

    // ASSERT_EQ(biconn.inverse_tree_map[1], 20);
    // ASSERT_EQ(biconn.tree_map[20].size(), 1);
    // ASSERT_EQ(biconn.tree_map[20][1], 0);
}

TEST(BiconnectedPartition, DoubleStarGraph) {
    NetworKit::Graph g = mcFileToGraph("test/data/double_star.mc");

    sms::BiconnectedPartition biconn(g);
    biconn.run();

    ASSERT_EQ(biconn.components().size(), 40);
    // ASSERT_EQ(biconn.bridges().size(), 40);

    // ASSERT_EQ(biconn.tree.numberOfNodes(), 61);
    // ASSERT_EQ(biconn.tree.numberOfEdges(), 60);
    // ASSERT_EQ(biconn.tree_map[5].size(), 2);

    // ASSERT_EQ(biconn.inverse_tree_map[1], 40);
    // ASSERT_EQ(biconn.tree_map[40].size(), 1);
    // ASSERT_EQ(biconn.tree_map[40][0], 1);

    // ASSERT_EQ(biconn.inverse_tree_map[1], 40);
    // ASSERT_EQ(biconn.tree_map[41].size(), 1);
    // ASSERT_EQ(biconn.tree_map[41][1], 0);
}