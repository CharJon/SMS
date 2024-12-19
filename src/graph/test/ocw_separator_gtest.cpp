#include <gtest/gtest.h>

#include "sms/graph/ocw_seperator.hpp"

TEST(OcwSeparator, SimpleNoPath) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 0.4);
    g.addEdge(1, 2, 0.4);
    g.addEdge(2, 0, 0.4);
    g.indexEdges();

    sms::OcwSeparator ocwSeparator(g);
    ocwSeparator.updateWeights(0, 1, 0.3);
    ocwSeparator.updateWeights(1, 2, 0.3);
    ocwSeparator.updateWeights(2, 0, 0.3);

    auto res = ocwSeparator.getViolatedOcws(0);
    EXPECT_EQ(res.size(), 0);
}

TEST(OcwSeparator, MostViolatedSimplePath) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 0.4);
    g.addEdge(1, 2, 0.4);
    g.addEdge(2, 0, 0.4);
    g.indexEdges();

    sms::OcwSeparator ocwSeparator(g);
    ocwSeparator.updateWeights(0, 1, 0.1);
    ocwSeparator.updateWeights(1, 2, 0.4);
    ocwSeparator.updateWeights(2, 0, 0.1);

    auto res = ocwSeparator.getMostViolatedOcw(0);
    ASSERT_TRUE(res.has_value());
    EXPECT_EQ(res.value().size(), 3);
}

TEST(OcwSeparator, SimplePath) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 0.4);
    g.addEdge(1, 2, 0.4);
    g.addEdge(2, 0, 0.4);
    g.indexEdges();

    sms::OcwSeparator ocwSeparator(g);
    ocwSeparator.updateWeights(0, 1, 0.9);
    ocwSeparator.updateWeights(1, 2, 0.9);
    ocwSeparator.updateWeights(2, 0, 0.9);

    auto res = ocwSeparator.getViolatedOcws(0);
    ASSERT_GE(res.size(), 1);
    ASSERT_TRUE(res[0].isValid(g));
}

TEST(OcwSeparator, ExtractSimple1) {
    NetworKit::Graph g(4);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(0, 3);
    g.indexEdges();

    sms::OcwSeparator ocws(g);
    ocws.updateWeights(0, 1, 1.);
    ocws.updateWeights(0, 2, 1.);
    ocws.updateWeights(1, 2, 1.);
    ocws.updateWeights(0, 3, 0.1);

    sms::OddClosedWalk ocw(3);
    ocw.addEdge(0, true);
    ocw.addEdge(1, false);
    ocw.addEdge(2, false);
    ocw.addEdge(0, false);
    ocw.addEdge(3, true);

    auto res = ocws.simpleOddCycles(ocw);
    EXPECT_NE(res.size(), 0);
}

TEST(OcwSeparator, Extract1) {
    NetworKit::Graph g(4);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.indexEdges();

    sms::OcwSeparator ocws(g);
    ocws.updateWeights(0, 1, 0.5);
    ocws.updateWeights(0, 2, 0.3);
    ocws.updateWeights(0, 3, 0.1);
    ocws.updateWeights(1, 2, 0.1);
    ocws.updateWeights(2, 3, 0.1);

    sms::OddClosedWalk ocw(0);
    ocw.addEdge(1, false);
    ocw.addEdge(2, true);
    ocw.addEdge(3, true);
    ocw.addEdge(0, true);

    ocws.simpleOddCycles(ocw);

    auto res = ocws.extractCycles(ocw);
    ASSERT_EQ(res.size(), 2);
    EXPECT_EQ(res[0].size(), 3);
    EXPECT_EQ(res[1].size(), 3);
}

TEST(OcwSeparator, Extract2) {
    NetworKit::Graph g(7);
    g.addEdge(0, 1);
    g.addEdge(0, 6);
    g.addEdge(1, 2);
    g.addEdge(1, 6);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(2, 5);
    g.addEdge(3, 4);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.indexEdges();

    sms::OcwSeparator ocws(g);
    ocws.updateWeights(0, 1, 0.1);
    ocws.updateWeights(0, 6, 0.1);
    ocws.updateWeights(1, 2, 0.9);
    ocws.updateWeights(1, 6, 0.5);
    ocws.updateWeights(2, 3, 0.9);
    ocws.updateWeights(2, 4, 0.1);
    ocws.updateWeights(2, 5, 0.1);
    ocws.updateWeights(3, 4, 0.1);
    ocws.updateWeights(4, 5, 0.1);
    ocws.updateWeights(5, 6, 0.9);

    sms::OddClosedWalk ocw(0);
    ocw.addEdge(1, true);
    ocw.addEdge(2, false);
    ocw.addEdge(3, false);
    ocw.addEdge(4, true);
    ocw.addEdge(5, true);
    ocw.addEdge(6, false);
    ocw.addEdge(0, true);

    ocws.simpleOddCycles(ocw);
    auto res = ocws.extractCycles(ocw);
    ASSERT_NE(res.size(), 0);
    EXPECT_EQ(res.size(), 3);
}
