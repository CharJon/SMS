#include <gtest/gtest.h>

#include "sms/instance/frustration_index.hpp"

TEST(FrustrationIndex, createClass1) {
    NetworKit::Graph g(0, true);
    sms::FrustrationIndex frustrationIndex(g);
    ASSERT_EQ(frustrationIndex.getNumberOfNodes(), 0);
    ASSERT_EQ(frustrationIndex.getNumberOfEdges(), 0);
}

TEST(FrustrationIndex, createClass2) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 0, -1.0);

    sms::FrustrationIndex frustrationIndex(g);
    ASSERT_EQ(frustrationIndex.getNumberOfNodes(), 3);
    ASSERT_EQ(frustrationIndex.getNumberOfEdges(), 3);
    ASSERT_EQ(frustrationIndex.getSolutionValue({1, 1, 1}), 1);
    ASSERT_EQ(frustrationIndex.getSolutionValue({1, 0, 1}), 3);
    ASSERT_EQ(frustrationIndex.getSolutionValue({0, 0, 0}), 1);

    ASSERT_TRUE(frustrationIndex.isConsistent());

    auto info = frustrationIndex.getInstanceInformation();

    ASSERT_EQ(info["num nodes"].get<int>(), 3);
    ASSERT_EQ(info["num edges"].get<int>(), 3);
    ASSERT_EQ(info["num negative edges"].get<int>(), 1);
    ASSERT_EQ(info["num positive edges"].get<int>(), 2);

    NetworKit::Graph g2(3, true);
    g2.addEdge(0, 1, 1.0);
    g2.addEdge(1, 2, 0.0);
    g2.addEdge(2, 0, -1.0);

    sms::FrustrationIndex frustrationIndex2(g2);

    ASSERT_FALSE(frustrationIndex2.isConsistent());
}

TEST(FrustrationIndex, createClass3) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 0, -2.0);

    sms::FrustrationIndex frustrationIndex(g);
    ASSERT_FALSE(frustrationIndex.isConsistent());
}
