#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/fast_vertex_separators.hpp"

TEST(FastVertexSeparator, Line) {
    auto g = NetworKit::Graph(5, true, false);
    g.addEdge(0, 1);
    g.addEdge(0, 4);
    g.addEdge(1, 2);
    g.addEdge(1, 3);

    auto sep = sms::FastVertexSeparator(g);
    sep.run(3);
    EXPECT_EQ(sep.separatedSets()[0].size(), 2);
}

TEST(FastVertexSeparator, CorrectNeighborhood) {
    auto g = NetworKit::Graph(10, true, false);
    g.addEdge(0, 1);
    g.addEdge(0, 4);
    g.addEdge(0, 3);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(1, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 4);
    g.addEdge(2, 5);
    g.addEdge(3, 4);
    g.addEdge(3, 6);
    g.addEdge(3, 7);
    g.addEdge(4, 5);
    g.addEdge(4, 6);
    g.addEdge(4, 7);
    g.addEdge(4, 8);
    g.addEdge(5, 7);
    g.addEdge(5, 8);
    g.addEdge(6, 7);
    g.addEdge(7, 8);

    auto sep = sms::FastVertexSeparator(g);
    sep.run(6);
    for (auto s : sep.separators()) {
        EXPECT_EQ(s.size(), 3);
    }
}
