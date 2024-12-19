#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs_statistics.hpp"

TEST(GraphStatistics, DegeneracySimple) {
    NetworKit::Graph g(6, false, false);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 5);
    g.addEdge(5, 1);
    g.addEdge(5, 2);
    g.addEdge(5, 3);

    EXPECT_EQ(sms::degeneracy(g, 10), 2);
}

TEST(GraphStatistics, DegeneracyTree) {
    NetworKit::Graph g(20, false, false);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 5);

    g.addEdge(6, 7);
    g.addEdge(6, 8);
    g.addEdge(8, 9);
    g.addEdge(8, 10);
    g.addEdge(10, 11);
    g.addEdge(7, 12);
    g.addEdge(7, 13);
    g.addEdge(7, 14);
    g.addEdge(7, 15);
    g.addEdge(7, 16);
    g.addEdge(16, 17);
    g.addEdge(17, 18);
    g.addEdge(17, 19);

    EXPECT_EQ(sms::degeneracy(g, 3), 1);
}

TEST(GraphStatistics, DegeneracyGrid) {
    NetworKit::Graph g(11, false, false);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(1, 2);
    g.addEdge(1, 4);
    g.addEdge(2, 3);
    g.addEdge(2, 5);

    g.addEdge(4, 7);
    g.addEdge(4, 5);
    g.addEdge(5, 8);
    g.addEdge(5, 6);

    g.addEdge(7, 8);
    g.addEdge(8, 9);
    g.addEdge(6, 9);
    g.addEdge(6, 3);

    g.addEdge(7, 10);
    g.addEdge(8, 10);
    g.addEdge(9, 10);

    EXPECT_EQ(sms::degeneracy(g, 10), 3);
    EXPECT_EQ(sms::degeneracy(g, 2), std::nullopt);
}

TEST(GraphStatistics, DegeneracyGrid2) {
    NetworKit::Graph g(11, false, false);
    for (int i = 1; i < 10; ++i) {
        g.addEdge(0, i);
        g.addEdge(10, i);
    }
    g.addEdge(1, 2);
    g.addEdge(1, 4);
    g.addEdge(2, 3);
    g.addEdge(2, 5);

    g.addEdge(4, 7);
    g.addEdge(4, 5);
    g.addEdge(5, 8);
    g.addEdge(5, 6);

    g.addEdge(7, 8);
    g.addEdge(8, 9);
    g.addEdge(6, 9);
    g.addEdge(6, 3);

    EXPECT_EQ(sms::degeneracy(g, 10), 4);
    EXPECT_EQ(sms::degeneracy(g, 2), std::nullopt);
}

TEST(GraphStatistics, DegeneracyClique) {
    NetworKit::Graph g(11, false, false);

    for (int i = 0; i < 11; ++i) {
        for (int j = i + 1; j < 11; ++j) {
            g.addEdge(i, j);
        }
    }

    EXPECT_EQ(sms::degeneracy(g, 10), 10);
    EXPECT_FALSE(sms::degeneracy(g, 2).has_value());
}

TEST(GraphStatistics, DegeneracyCliqueWithAppendix) {
    NetworKit::Graph g(14, false, false);

    for (int i = 0; i < 11; ++i) {
        for (int j = i + 1; j < 11; ++j) {
            g.addEdge(i, j);
        }
    }
    g.addEdge(10, 11);
    g.addEdge(11, 12);
    g.addEdge(12, 13);

    EXPECT_EQ(sms::degeneracy(g, 15), 10);
    EXPECT_FALSE(sms::degeneracy(g, 2).has_value());
}

TEST(GraphStatistics, DegeneracyCliqueWithClique1) {
    NetworKit::Graph g(25, false, false);

    for (int i = 0; i < 11; ++i) {
        for (int j = i + 1; j < 11; ++j) {
            g.addEdge(i, j);
        }
    }

    for (int i = 11; i < 25; ++i) {
        for (int j = i + 1; j < 25; ++j) {
            g.addEdge(i, j);
        }
    }

    EXPECT_EQ(sms::degeneracy(g, 25), 13);
    EXPECT_FALSE(sms::degeneracy(g, 10).has_value());
}

TEST(GraphStatistics, DegeneracyCliqueWithClique2) {
    NetworKit::Graph g(25, false, false);

    for (int i = 0; i < 11; ++i) {
        for (int j = i + 1; j < 11; ++j) {
            g.addEdge(i, j);
        }
    }
    g.addEdge(10, 11);
    for (int i = 11; i < 25; ++i) {
        for (int j = i + 1; j < 25; ++j) {
            g.addEdge(i, j);
        }
    }

    EXPECT_EQ(sms::degeneracy(g, 25), 13);
    EXPECT_FALSE(sms::degeneracy(g, 10).has_value());
}
