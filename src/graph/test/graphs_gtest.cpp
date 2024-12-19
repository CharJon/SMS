#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "networkit/io/EdgeListReader.hpp"

#include "sms/auxiliary/math.hpp"
#include "sms/graph/graphs.hpp"
#include "sms/graph/graphs_statistics.hpp"
#include "sms/io/io.hpp"

void test_subgraph(NetworKit::Graph &g, const std::vector<NetworKit::node> &nodes) {
    auto [inducedSubgraph, sub2orig, orig2sub] = sms::inducedSubgraphCompact(g, nodes);

    unsigned int c = 0;
    for (auto i : nodes) {
        ASSERT_EQ(orig2sub[i], c);
        c++;
        ASSERT_EQ(i, sub2orig[orig2sub[i]]);
    }

    for (unsigned int i = 0; i < nodes.size(); i++) {
        ASSERT_EQ(sub2orig[i], nodes[i]);
    }

    for (auto e : inducedSubgraph.edgeWeightRange())
        ASSERT_EQ(e.weight, g.weight(sub2orig[e.u], sub2orig[e.v]));

    auto nkSub = NetworKit::GraphTools::subgraphFromNodes(g, nodes.begin(), nodes.end(), true);

    ASSERT_EQ(nkSub.numberOfNodes(), inducedSubgraph.numberOfNodes());
    ASSERT_EQ(nkSub.numberOfEdges(), inducedSubgraph.numberOfEdges());
    ASSERT_EQ(nkSub.totalEdgeWeight(), inducedSubgraph.totalEdgeWeight());
}

TEST(AuxGraphs, InducedSubgraphSimpleTest) {
    NetworKit::Graph g(4);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 3);
    g.addEdge(2, 3);

    std::vector<NetworKit::node> nodes = {1, 3};
    auto [s, _, __] = sms::inducedSubgraphCompact(g, nodes);
    EXPECT_EQ(s.numberOfNodes(), 2);
    EXPECT_EQ(s.numberOfEdges(), 1);
    EXPECT_EQ(s.upperNodeIdBound(), 2);
    EXPECT_TRUE(s.hasEdge(0, 1));
}

TEST(AuxGraphs, InducedSubgraphComplexTest) {
    NetworKit::Graph g = mcFileToGraph("test/data/ego-facebook.mc");

    NetworKit::count subgraphSize = g.numberOfNodes() * 2 / 3;
    std::vector<NetworKit::node> nodes;

    for (const auto u : g.nodeRange()) {
        nodes.emplace_back(u);
        if (nodes.size() >= subgraphSize)
            break;
    }

    test_subgraph(g, nodes);
}

TEST(AuxGraphs, InducedSubgraphComplexTest2) {
    NetworKit::Graph g = mcFileToGraph("test/data/subgraph_test.mc");

    std::vector<NetworKit::node> nodes = {1, 2, 3};

    test_subgraph(g, nodes);
}

TEST(AuxGraphs, InducedSubgraphComplexTest3) {
    NetworKit::EdgeListReader el(' ', 0);
    NetworKit::Graph g = el.read("test/data/ia-wiki-Talk.el");

    NetworKit::count subgraphSize = g.numberOfNodes() * 2 / 3;
    std::vector<NetworKit::node> nodes(subgraphSize);

    for (unsigned int i = 0; i < subgraphSize; i++) {
        nodes[i] = i;
    }

    test_subgraph(g, nodes);
}

TEST(AuxGraphs, InducedSubgraphComplexTest3Nk) {
    NetworKit::EdgeListReader el(' ', 0);
    NetworKit::Graph g = el.read("test/data/ia-wiki-Talk.el");

    NetworKit::count subgraphSize = g.numberOfNodes() * 2 / 3;
    std::vector<NetworKit::node> nodes(subgraphSize);

    for (unsigned int i = 0; i < subgraphSize; i++) {
        nodes[i] = i;
    }

    test_subgraph(g, nodes);
}

TEST(AuxGraphs, compactGraphShuffled1) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 0.7);
    g.addEdge(1, 3, 0.7);
    g.removeNode(2);

    auto compactG = sms::compactGraph(g, 4);

    EXPECT_EQ(compactG.compactGraph.numberOfSelfLoops(), g.numberOfSelfLoops());
    EXPECT_EQ(compactG.compactGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(compactG.compactGraph.numberOfEdges(), g.numberOfEdges());
}

TEST(AuxGraphs, complementGraphLine) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 3, 1);
    ASSERT_EQ(g.numberOfNodes(), 4);

    auto complementGraph = sms::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 3);
}

TEST(AuxGraphs, complementGraphCycle1) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 4, 1);
    g.addEdge(4, 3, 1);
    g.addEdge(3, 0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = sms::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 5);
}

TEST(AuxGraphs, complementGraphCycle2) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 4, 1);
    g.addEdge(4, 1, -1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(3, 0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = sms::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 5);
}

TEST(AuxGraphs, complementGraphClique) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v);
        }
    }
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = sms::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 0);
}

TEST(AuxGraphs, complementGraphNearClique) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v);
        }
    }
    g.removeEdge(0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = sms::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 1);
}

TEST(AuxGraphs, neighborhoodAlphaConnected) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 47);
    g.addEdge(0, 2, 2);
    g.addEdge(1, 2, 4);
    g.addEdge(0, 3, -5);
    g.addEdge(1, 3, -10);
    g.addEdge(0, 4, 7);
    g.addEdge(1, 4, 14);
    ASSERT_EQ(g.numberOfNodes(), 5);

    std::vector<NetworKit::edgeweight> aux(g.upperNodeIdBound(), 0.0);
    EXPECT_EQ(sms::neighborhoodAlpha(g, 0, 1, true, aux), 0.5);
    EXPECT_EQ(sms::neighborhoodAlpha(g, 0, 1, false, aux), 0);
}

TEST(AuxGraphs, neighborhoodAlphaDisconnected) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 2, 2);
    g.addEdge(1, 2, -3);
    g.addEdge(0, 3, 3);
    g.addEdge(1, 3, -4.5);
    g.addEdge(0, 4, -5);
    g.addEdge(1, 4, 7.5);
    ASSERT_EQ(g.numberOfNodes(), 5);

    std::vector<NetworKit::edgeweight> aux(g.upperNodeIdBound(), 0.0);
    EXPECT_EQ(sms::neighborhoodAlpha(g, 0, 1, true, aux), -2. / 3.);
    EXPECT_EQ(sms::neighborhoodAlpha(g, 0, 1, false, aux), -2. / 3.);
}

TEST(AuxGraphs, edgeWeightDivisorInt) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v, 2.);
        }
    }
    g.removeEdge(0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto commonDivisor = sms::edgeWeightDivisor(g);
    EXPECT_EQ(commonDivisor, 2.0);
}

TEST(AuxGraphs, edgeWeightDivisorFloatMinMax) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v, 1.5);
        }
    }
    g.setWeight(0, 1, -1.5);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto commonDivisor = sms::edgeWeightDivisor(g);
    EXPECT_EQ(commonDivisor, 1.5);
}

TEST(AuxGraphs, edgeWeightDivisorFloat) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v, 1.7);
        }
    }
    g.setWeight(0, 1, -1.71);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto commonDivisor = sms::edgeWeightDivisor(g);
    EXPECT_EQ(commonDivisor, 1.);
}

TEST(AuxGraphs, edgeWeightDivisorDegree) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v, 1.7);
        }
    }
    g.setWeight(0, 1, -1.71);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto divisor = sms::degreeBasedScaling(g);
    EXPECT_EQ(divisor, 1.);
}

TEST(AuxGraphs, edgeWeightDivisorEvenDegreeFloat) {
    NetworKit::Graph g(8, true, false);
    for (NetworKit::node u = 0; u < 4; u++) {
        for (NetworKit::node v = 4; v < 8; v++) {
            g.addEdge(u, v, 1.7);
        }
    }
    ASSERT_EQ(g.numberOfNodes(), 8);
    ASSERT_EQ(g.numberOfEdges(), 16);

    auto divisor = sms::degreeBasedScaling(g);
    EXPECT_EQ(divisor, 1.);
}

TEST(AuxGraphs, edgeWeightDivisorEvenDegreeInt) {
    NetworKit::Graph g(8, true, false);
    for (NetworKit::node u = 0; u < 4; u++) {
        for (NetworKit::node v = 4; v < 8; v++) {
            g.addEdge(u, v, 3.0);
        }
    }
    g.removeEdge(0, 4);
    g.removeEdge(0, 5);
    g.removeEdge(1, 4);
    g.removeEdge(1, 5);
    g.setWeight(0, 6, 5.0);

    ASSERT_EQ(g.numberOfNodes(), 8);
    ASSERT_EQ(g.numberOfEdges(), 12);

    auto divisor = sms::degreeBasedScaling(g);
    EXPECT_EQ(divisor, 2.);
}

TEST(AuxGraphs, edgeWeightDivisorEvenAndOddDegreeInt) {
    std::cout << int(-3.0) % 2 << '\n';
    NetworKit::Graph g(9, true, false);
    for (NetworKit::node u = 0; u < 4; u++) {
        for (NetworKit::node v = 4; v < 8; v++) {
            g.addEdge(u, v, 3.0);
        }
    }
    g.removeEdge(0, 4);
    g.removeEdge(0, 5);
    g.removeEdge(1, 4);
    g.removeEdge(1, 5);
    g.setWeight(0, 6, 5.0);
    g.addEdge(0, 8, 2.0);

    ASSERT_EQ(g.numberOfNodes(), 9);
    ASSERT_EQ(g.numberOfEdges(), 13);

    auto divisor = sms::degreeBasedScaling(g);
    EXPECT_EQ(divisor, 2.);
}

TEST(Graph, BasicTest1) {
    NetworKit::Graph graph(3, true);

    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 1);
    graph.addEdge(0, 2, 0.5);

    EXPECT_FALSE(sms::integerWeightedDegree(graph));
}

TEST(Graph, BasicTest2) {
    NetworKit::Graph graph(3, true);

    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 1);
    graph.addEdge(0, 2, 1);

    EXPECT_TRUE(sms::integerWeightedDegree(graph));
}

TEST(Graph, BasicTest3) {
    NetworKit::Graph graph(3, true);

    graph.addEdge(0, 1, 0.5);
    graph.addEdge(1, 2, 0.5);
    graph.addEdge(0, 2, 0.5);

    EXPECT_TRUE(sms::integerWeightedDegree(graph));
}

TEST(Graph, ScalingSucess) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 2, 2.4);
    g.addEdge(1, 2, -3.2);
    g.addEdge(0, 3, 3.2);
    g.addEdge(1, 3, -4);
    g.addEdge(0, 4, -5.4);
    g.addEdge(1, 4, 7.8);

    auto scalar = sms::integerScalar(g);
    auto graph = sms::scaleToInt(g, scalar);

    for (auto e : graph.edgeWeightRange()) {
        EXPECT_TRUE(sms::isInteger(e.weight));
    }
}

TEST(Graph, MinWeightScale) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 2, 2.4);
    g.addEdge(1, 2, -3.2);
    g.addEdge(0, 3, 3.2);
    g.addEdge(1, 3, -4);
    g.addEdge(0, 4, -5.4);
    g.addEdge(1, 4, 1. / 500);

    auto scalar = sms::integerScalar(g);

    EXPECT_EQ(scalar, 500);

    auto graph = sms::scaleToInt(g, scalar);

    for (auto e : graph.edgeWeightRange()) {
        EXPECT_TRUE(sms::isInteger(e.weight));
    }
}

TEST(Graph, ScalingFail) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 2, 2.40385500275729);
    g.addEdge(1, 2, -3.324085202);
    g.addEdge(0, 3, 3.258293);
    g.addEdge(1, 3, -4.510482);
    g.addEdge(0, 4, -5.074047);
    g.addEdge(1, 4, 7.550739540);

    auto scalar1 = sms::integerScalar(g, 10e5);
    EXPECT_EQ(scalar1, 0);

    auto scalar2 = sms::integerScalar(g);
    auto graph = sms::scaleToInt(g, scalar2);

    for (auto e : graph.edgeWeightRange()) {
        EXPECT_TRUE(sms::isInteger(e.weight));
    }
}