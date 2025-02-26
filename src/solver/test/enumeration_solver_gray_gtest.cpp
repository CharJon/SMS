#include <vector>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"
#include "networkit/io/EdgeListReader.hpp"

#include "sms/io/io.hpp"
#include "sms/solver/enumeration_solver_gray.hpp"

NetworKit::Graph squareGraph() {
    auto er = NetworKit::EdgeListReader();
    auto g = er.read("test/data/square.wel");
    assert(g.isWeighted());
    return g;
}

NetworKit::Graph cliqueGraph(int n) {
    assert(n > 0);
    auto g = NetworKit::Graph(n, true);

    for (auto u : g.nodeRange()) {
        for (auto v : g.nodeRange()) {
            if (u < v)
                g.addEdge(u, v, 1);
        }
    }

    return g;
}

TEST(Enumeration, SolvingSimple) {
    NetworKit::Graph g = mcFileToGraph("test/data/square.mc");
    g.restoreNode(0);

    sms::EnumerationSolverGray solver(g);

    solver.run();

    ASSERT_NEAR(solver.bestSolutionValue(), 1.1, 1e-6);
}

TEST(Enumeration, SolvingComplex) {
    NetworKit::Graph g = mcFileToGraph("test/data/easy_20_vertices.mc");
    g.restoreNode(0);
    sms::EnumerationSolverGray solver(g);
    solver.run();
    ASSERT_EQ(solver.bestSolutionValue(), 59);
}

class EnumerationParameterizedTestFixture : public ::testing::TestWithParam<int> {};

TEST_P(EnumerationParameterizedTestFixture, SolvingCliques) {
    const int dim = GetParam();
    auto g = cliqueGraph(dim);
    sms::EnumerationSolverGray solver(g);
    solver.run();

    ASSERT_EQ(solver.bestSolutionValue(), (dim / 2) * ((dim + 1) / 2));
}

INSTANTIATE_TEST_SUITE_P(AllSizes, EnumerationParameterizedTestFixture, ::testing::Values(2, 5, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23));


