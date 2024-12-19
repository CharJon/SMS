#include <vector>

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/small_ccs.hpp"
#include "sms/io/io.hpp"

TEST(smallCCs, sortFourCycles) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 0);

    g.addEdge(3, 4);

    auto r = sms::sortFourCycle({0, 1, 2, 3}, g);
    sms::fourcycle c = {1, 0, 3, 2};
    ASSERT_EQ(r, c);

    NetworKit::Graph g_2(4);
    g_2.addEdge(3, 1);
    g_2.addEdge(1, 0);
    g_2.addEdge(2, 0);
    g_2.addEdge(3, 2);

    r = sms::sortFourCycle({0, 1, 2, 3}, g_2);
    c = {1, 0, 2, 3};
    ASSERT_EQ(r, c);
}

TEST(smallCCs, sortTriangles) {
    auto r = sms::sortTriangle({1, 2, 3});
    sms::triangle c = {2, 1, 3};
    ASSERT_EQ(r, c);

    r = sms::sortTriangle({8, 9, 3});
    c = {8, 3, 9};
    ASSERT_EQ(r, c);
}

TEST(smallCCs, trippleSquareGraph) {
    NetworKit::Graph g = mcFileToGraph("test/data/tripleSquareTriangle.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), 2);
    ASSERT_EQ(ccs.fourHoles.size(), 2);

    EXPECT_THAT(ccs.triangles, ::testing::Contains(::testing::ElementsAre(6, 5, 7)));
    EXPECT_THAT(ccs.triangles, ::testing::Contains(::testing::ElementsAre(7, 5, 8)));

    EXPECT_THAT(ccs.fourHoles, ::testing::Contains(::testing::ElementsAre(2, 1, 4, 3)));
    EXPECT_THAT(ccs.fourHoles, ::testing::Contains(::testing::ElementsAre(4, 3, 6, 5)));
}

TEST(smallCCs, noTriangle) {
    NetworKit::Graph g = mcFileToGraph("test/data/square.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), 0);
    ASSERT_EQ(ccs.fourHoles.size(), 1);
}

TEST(smallCCs, noFourHoles) {
    NetworKit::Graph g = mcFileToGraph("test/data/two_triangles.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), 2);
    ASSERT_EQ(ccs.fourHoles.size(), 0);
}

TEST(smallCCs, nothing) {
    NetworKit::Graph g = mcFileToGraph("test/data/double_star.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), 0);
    ASSERT_EQ(ccs.fourHoles.size(), 0);
}

TEST(smallCCs, caNetsceince) {
    NetworKit::Graph g = mcFileToGraph("test/data/ca-netscience.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());

    ASSERT_EQ(ccs.triangles.size(), 921);
}

TEST(smallCCs, egoFacebook) {
    NetworKit::Graph g = mcFileToGraph("test/data/ego-facebook.mc");
    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), 273 / 3);
}

TEST(smallCCs, complete) {
    int c = 100;
    NetworKit::Graph g(c);

    for (int i = 0; i < c; i++) {
        for (int j = i + 1; j < c; j++) {
            g.addEdge(i, j);
        }
    }

    g.indexEdges();

    sms::SmallChordlessCycles ccs(g);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());
    ASSERT_EQ(ccs.triangles.size(), c * (c - 1) * (c - 2) / 6);
}