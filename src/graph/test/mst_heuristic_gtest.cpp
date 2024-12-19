#include <gtest/gtest.h>

#include "sms/graph/mst_heuristic.hpp"
#include "sms/instance/mc_solution.hpp"

TEST(MSTHeuristic, PrimalHeurPath) {
    NetworKit::Graph g(4, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.5);
    heur.updateWeights(1, 2, 0.5);
    heur.updateWeights(2, 3, 0.5);

    heur.computeSpanningTree();

    ASSERT_EQ(sms::solutionValue(g, heur.getPrimalSolution(0)), 0);
}

TEST(MSTHeuristic, PrimalHeurSimple) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(1, 4, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.7);
    heur.updateWeights(1, 2, 0.4);
    heur.updateWeights(1, 3, 0.6);
    heur.updateWeights(3, 4, 0.49);
    heur.updateWeights(1, 4, 0.7);

    heur.computeSpanningTree();

    ASSERT_EQ(sms::solutionValue(g, heur.getPrimalSolution(0)), 3);

    // auto ocs = heur.getViolatedOCs();
    // ASSERT_EQ(ocs.size(), 0);
}

TEST(MSTHeuristic, ConsSeparatorSimple1) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(1, 4, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.7);
    heur.updateWeights(1, 2, 0.4);
    heur.updateWeights(1, 3, 0.6);
    heur.updateWeights(3, 4, 0.8);
    heur.updateWeights(1, 4, 0.7);

    heur.computeSpanningTree();
    auto path = heur.getSTPath(1, 3);

    ASSERT_EQ(path.size(), 3);

    // auto ocs = heur.getViolatedOCs();
    // ASSERT_EQ(ocs.size(), 1);
    // ASSERT_EQ(ocs[0].size(), 3);
    // ASSERT_TRUE(ocs[0].isValid(g));
}

TEST(MSTHeuristic, ConsSeparatorTriangleRemove) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(0, 2, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.9);
    heur.updateWeights(1, 2, 0.6);
    heur.updateWeights(0, 2, 0.1);

    heur.computeSpanningTree();

    // auto ocs = heur.getViolatedOCs();
    // ASSERT_EQ(ocs.size(), 1);
    // ASSERT_TRUE(ocs[0].isValid(g));
    // ASSERT_EQ(ocs[0].size(), 3);
}

TEST(MSTHeuristic, ConsSeparatorTriangleAdd) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(0, 2, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.9);
    heur.updateWeights(1, 2, 0.9);
    heur.updateWeights(0, 2, 0.4);

    heur.computeSpanningTree();
    auto path = heur.getSTPath(0, 2);

    ASSERT_EQ(path.size(), 3);

    // auto ocs = heur.getViolatedOCs();
    // ASSERT_EQ(ocs.size(), 1);
    // ASSERT_TRUE(ocs[0].isValid(g));
    // ASSERT_EQ(ocs[0].size(), 3);
}

TEST(MSTHeuristic, ConsSeparatorComplex) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(1, 4, 1);
    g.addEdge(2, 3, 1);

    sms::MSTHeuristic heur(&g);

    heur.updateWeights(0, 1, 0.7);
    heur.updateWeights(1, 2, 0.39);
    heur.updateWeights(1, 3, 0.6);
    heur.updateWeights(3, 4, 0.8);
    heur.updateWeights(1, 4, 0.7);
    heur.updateWeights(2, 3, 0.9);

    heur.computeSpanningTree();

    // auto ocs = heur.getViolatedOCs();
    //  ASSERT_EQ(ocs.size(), 2);
}
