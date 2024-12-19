#include <chrono>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/two_separator.hpp"
#include "sms/io/io.hpp"

TEST(TwoSeparator, BasicSquareTest) {
    NetworKit::Graph g(4);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 0);

    sms::TwoSeparator sep(g);

    sep.run();

    ASSERT_TRUE(sep.hasRun());

    ASSERT_EQ(sep.separators.size(), 2);

    ASSERT_EQ(sep.separators[0], std::make_tuple(0, 2));
    ASSERT_EQ(sep.separators[1], std::make_tuple(1, 3));
}

TEST(TwoSeparator, MediumSizeTest) {
    NetworKit::Graph g(10);
    for (unsigned int v = 0; v < 5; v++) {
        for (unsigned int u = 0; u < 5; u++) {
            if (u < v) {
                g.addEdge(v, u);
                g.addEdge(v + 5, u + 5);
            }
        }
    }

    g.addEdge(0, 6);
    g.addEdge(1, 5);

    ASSERT_EQ(g.numberOfEdges(), 22);

    sms::TwoSeparator sep(g);

    sep.run();

    ASSERT_TRUE(sep.hasRun());

    ASSERT_EQ(sep.separators.size(), 4);

    ASSERT_EQ(sep.separators[0], std::make_tuple(0, 1));
    ASSERT_EQ(sep.separators[1], std::make_tuple(0, 5));
    ASSERT_EQ(sep.separators[2], std::make_tuple(1, 6));
    ASSERT_EQ(sep.separators[3], std::make_tuple(5, 6));
}

TEST(TwoSeparator, ThreeComponents) {
    NetworKit::Graph g(17);
    for (unsigned int v = 0; v < 5; v++) {
        for (unsigned int u = 0; u < 5; u++) {
            if (u < v) {
                g.addEdge(v, u);
                g.addEdge(v + 5, u + 5);
                g.addEdge(v + 10, u + 10);
            }
        }
    }

    g.addEdge(15, 0);
    g.addEdge(15, 5);
    g.addEdge(15, 10);
    g.addEdge(15, 1);
    g.addEdge(15, 6);
    g.addEdge(15, 11);

    g.addEdge(16, 0);
    g.addEdge(16, 5);
    g.addEdge(16, 10);
    g.addEdge(16, 4);
    g.addEdge(16, 9);
    g.addEdge(16, 14);

    ASSERT_EQ(g.numberOfEdges(), 42);

    sms::TwoSeparator sep(g);

    sep.run();

    ASSERT_TRUE(sep.hasRun());

    ASSERT_EQ(sep.separators.size(), 1);

    ASSERT_EQ(sep.separators[0], std::make_tuple(15, 16));
}

TEST(TwoSeparator, SimpleRunTimeTest) {
    // runtime on 22.05.2023 (before fixing use of DynConnectedComp): ~ 10 sec (coan-22)
    // runtime on 23.05.2023 (after fixing use of DynConnectedComp): ~ 0.6 sec (coan-22)

    NetworKit::Graph g(100);
    for (unsigned int v = 0; v < 50; v++) {
        for (unsigned int u = 0; u < 50; u++) {
            if (u < v) {
                g.addEdge(v, u);
                g.addEdge(v + 50, u + 50);
            }
        }
    }

    g.addEdge(0, 51);
    g.addEdge(1, 50);

    sms::TwoSeparator sep(g);

    auto start = std::chrono::high_resolution_clock::now();
    sep.run();
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    ASSERT_TRUE(sep.hasRun());

    ASSERT_EQ(sep.separators.size(), 4);
}