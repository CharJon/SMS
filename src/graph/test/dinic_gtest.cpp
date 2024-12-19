#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/dinic.hpp"

TEST(Dinic, Line) {
    auto g = NetworKit::Graph(3, true, false, true);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 2.0);

    auto dinic = sms::Dinic(g);
    dinic.run(0, 2);

    ASSERT_EQ(dinic.getFlowValue(), 1);
}

TEST(Dinic, Triangle) {
    auto g = NetworKit::Graph(3, true, false, true);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 2.0);
    g.addEdge(0, 2, 2.0);

    auto dinic = sms::Dinic(g);
    dinic.run(0, 2);

    ASSERT_EQ(dinic.getFlowValue(), 3);
}

TEST(Dinic, PentagonGraph) {
    auto g = NetworKit::Graph(5, true, false, true);
    g.addEdge(0, 1, 2.0);
    g.addEdge(0, 2, 5.0);
    g.addEdge(1, 2, 4.0);
    g.addEdge(1, 3, 7.0);
    g.addEdge(2, 3, 8.0);
    g.addEdge(2, 4, 2.0);
    g.addEdge(3, 4, 12.0);

    auto dinic = sms::Dinic(g);
    dinic.run(0, 4);

    ASSERT_EQ(dinic.getFlowValue(), 7);
}

TEST(Dinic, SendFlowBack) {
    auto g = NetworKit::Graph(7, true, false, true);
    g.addEdge(0, 1, 5.0);
    g.addEdge(0, 2, 8.0);

    g.addEdge(1, 3, 5.0);
    g.addEdge(2, 3, 8.0);
    g.addEdge(3, 4, 7.0);
    g.addEdge(6, 4, 5.0);
    g.addEdge(1, 5, 5.0);
    g.addEdge(5, 6, 5.0);

    auto dinic = sms::Dinic(g);
    dinic.run(0, 4);

    ASSERT_EQ(dinic.getFlowValue(), 7 + 5);

    auto sourceSet = dinic.getSourceSet();

    ASSERT_EQ(sourceSet.size(), 4);
}