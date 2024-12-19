#include <gtest/gtest.h>

#include "sms/instance/maxcut.hpp"

TEST(MaxCut, createClass1) {
    NetworKit::Graph g(0);
    sms::MaxCut maxcut(g);
    ASSERT_EQ(maxcut.getNumberOfVertices(), 0);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 0);
    ASSERT_EQ(maxcut.getScalingFactor(), 1);
    ASSERT_EQ(maxcut.getOffset(), 0);
}

TEST(MaxCut, createClass2) {
    NetworKit::Graph g(3);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);

    sms::MaxCut maxcut(g);
    ASSERT_EQ(maxcut.getNumberOfVertices(), 3);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 3);
    ASSERT_EQ(maxcut.getScalingFactor(), 1);
    ASSERT_EQ(maxcut.getOffset(), 0);

    auto info = maxcut.getInstanceInformation();

    ASSERT_EQ(info["num nodes"].get<int>(), 3);
    ASSERT_EQ(info["num edges"].get<int>(), 3);
    ASSERT_EQ(info["scaling factor"].get<double>(), 1);
    ASSERT_EQ(info["offset"].get<double>(), 0);
}
