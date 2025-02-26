#include <vector>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"
#include "networkit/io/EdgeListReader.hpp"

#include "sms/io/io.hpp"
#include "sms/solver/enumeration_solver.hpp"

NetworKit::Graph squareGraph() {
    auto er = NetworKit::EdgeListReader();
    auto g = er.read("test/data/square.wel");
    assert(g.isWeighted());
    return g;
}

NetworKit::Graph cliqueGraph(int n) {
    assert(n > 0);
    auto g = NetworKit::Graph(n);

    for (auto u : g.nodeRange()) {
        for (auto v : g.nodeRange()) {
            if (u < v)
                g.addEdge(u, v);
        }
    }

    return g;
}

TEST(Enumeration, CutValue) {
    NetworKit::Graph g = squareGraph();

    sms::EnumerationSolver solver(g);

    unsigned int cut = 0;
    ASSERT_EQ(solver.getCutValue(cut), 0);

    cut = 15;
    ASSERT_EQ(solver.getCutValue(cut), 0);

    cut = 1;
    ASSERT_EQ(solver.getCutValue(cut), -0.7);

    cut = 14;
    ASSERT_EQ(solver.getCutValue(cut), -0.7);

    cut = 9;
    ASSERT_EQ(solver.getCutValue(cut), 1.1);
}

TEST(Enumeration, SolvingSimple) {
    NetworKit::Graph g = mcFileToGraph("test/data/square.mc");
    g.restoreNode(0);

    sms::EnumerationSolver solver(g);

    solver.run();

    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingComplex) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolver solver(g);
    solver.run();
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}

class EnumerationParameterizedTestFixture : public ::testing::TestWithParam<int> {
};

TEST_P(EnumerationParameterizedTestFixture, SolvingCliques) {
    const int dim = GetParam();
    auto g = cliqueGraph(dim);
    sms::EnumerationSolver solver(g);
    solver.run();

    ASSERT_EQ(solver.bestSolutionValue(), (dim / 2) * ((dim + 1) / 2));
}

INSTANTIATE_TEST_SUITE_P(AllSizes, EnumerationParameterizedTestFixture,
                         ::testing::Values(2, 5, 10, 15, 16, 17, 18, 19, 20, 21));

TEST(Enumeration, SolvingComplexParallel) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolver solver(g);
    solver.run();
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}

TEST(Enumeration, SolvingSimpleFixedLast1) {
    NetworKit::Graph g = squareGraph();

    sms::EnumerationSolver solver(g);
    solver.run({3}, {true});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedLast2) {
    NetworKit::Graph g = squareGraph();

    sms::EnumerationSolver solver(g);
    solver.run({3}, {false});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedFirst1) {
    NetworKit::Graph g = squareGraph();

    sms::EnumerationSolver solver(g);
    solver.run({0}, {true});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedFirst2) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);

    solver.run({0}, {false});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedMulti1) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);

    solver.run({1, 2}, {false, true});
    ASSERT_EQ(solver.bestSolutionValue(), -1.4);
}

TEST(Enumeration, SolvingSimpleFixedMulti2) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);

    solver.run({1, 2}, {true, true});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedMulti3) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);
    solver.run({1, 2}, {false, false});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedAlmostAll1) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);

    solver.run({0, 2, 3}, {true, false, true});
    ASSERT_EQ(solver.bestSolutionValue(), 1.1);
}

TEST(Enumeration, SolvingSimpleFixedAlmostAll2) {
    NetworKit::Graph g = squareGraph();
    sms::EnumerationSolver solver(g);

    solver.run({0, 1, 2}, {true, true, true});
    ASSERT_EQ(solver.bestSolutionValue(), 0);
}

TEST(Enumeration, SolvingComplexFixed1) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolver solver(g);
    solver.run({1}, {true});
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}

TEST(Enumeration, SolvingComplexFixed2) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolver solver(g);
    solver.run({0, 2}, {true, false});
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}

TEST(Enumeration, SolvingComplexFixed3) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolver solver(g);
    solver.run({0, 2, 3}, {true, false, true});
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}