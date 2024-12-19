#include <vector>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/solver/heuristic_solver.hpp"

NetworKit::Graph cliqueGraph(int n) {
    assert(n > 0);
    auto g = NetworKit::Graph(n, true, false);

    for (auto u : g.nodeRange()) {
        for (auto v : g.nodeRange()) {
            if (u < v)
                g.addEdge(u, v);
        }
    }

    return g;
}

TEST(HeuristicSolver, CutValue) {
    NetworKit::Graph g = cliqueGraph(50);

    sms::HeuristicSolver solver(g);
    solver.run();
    EXPECT_GE(solver.bestSolutionValue(), 1);
}
