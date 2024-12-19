#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/instance/maxcut.hpp"
#include "sms/solver/gurobi_quadratic_solver.hpp"

TEST(GurobiQuadraticSolver, SimpleTest) {

    NetworKit::Graph g(4, true);

    g.addEdge(0, 1, -1.0);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 3, -1.0);
    g.addEdge(0, 3, 2.0);

    auto solver = sms::GurobiQuadraticSolver(g);
    solver.setSeed(42);
    solver.setThreads(1);
    solver.setTimelimit(10);
    solver.run();
    auto solution =  solver.bestSolutionValue();

    ASSERT_EQ(solution, 1);
}
