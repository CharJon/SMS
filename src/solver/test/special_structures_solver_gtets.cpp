#include <vector>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"
#include "sms/solver/special_structures_solver.hpp"

TEST(SpecialStructuresSolver, tinyTestKhalfhalfUnweighted) {
    NetworKit::Graph g(5);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(0, 4);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(0, 1);
    g.addEdge(2, 3);

    NetworKit::Graph gCop = g;

    sms::SpecialStructuresSolver solver(g);
    EXPECT_TRUE(solver.applicable());
    solver.run();

    EXPECT_EQ(solver.bestSolutionValue(), 6);

    auto bestSolution = solver.getBestSolution();

    int numEdgesInCut = 0;

    for (int i = 0; i < 5; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (int j = i + 1; j < 5; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(gCop.hasEdge(i, j));
                numEdgesInCut++;
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, 6);
}

TEST(SpecialStructuresSolver, tinyTestKhalfhalfWeighted) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 2, 2);
    g.addEdge(0, 3, 2);
    g.addEdge(0, 4, 2);
    g.addEdge(1, 2, 2);
    g.addEdge(1, 3, 2);
    g.addEdge(1, 4, 2);
    g.addEdge(0, 1, -2);
    g.addEdge(2, 3, 2);
    g.addEdge(2, 4, -2);
    g.addEdge(3, 4, -2);

    NetworKit::Graph gCop = g;

    sms::SpecialStructuresSolver solver(g);
    EXPECT_TRUE(solver.applicable());
    solver.run();

    EXPECT_EQ(solver.bestSolutionValue(), 12);

    auto bestSolution = solver.getBestSolution();

    int numEdgesInCut = 0;

    for (int i = 0; i < 5; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (int j = i + 1; j < 5; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(gCop.hasEdge(i, j));
                numEdgesInCut++;
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, 6);
}

TEST(SpecialStructuresSolver, notApplicableUnweighted) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);

    sms::SpecialStructuresSolver solver(g);
    EXPECT_FALSE(solver.applicable());
    EXPECT_ANY_THROW(solver.run());
}

TEST(SpecialStructuresSolver, largeApplicableKhalfhalfUnweighted) {
    uint n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<uint> diss(0, n);

    NetworKit::Graph g(n);

    int numMissingEdges = 0;
    for (uint i = 0; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            if (i < n / 2 && j >= n / 2) {
                g.addEdge(i, j);
            } else {
                if (diss(rng) % 2 == 0)
                    g.addEdge(i, j);
                else
                    numMissingEdges++;
            }
        }
    }

    EXPECT_EQ(numMissingEdges, sms::unweightedComplementGraph(g).numberOfEdges());

    sms::SpecialStructuresSolver solver(g);
    EXPECT_TRUE(solver.applicable());
    solver.run();

    EXPECT_EQ(solver.bestSolutionValue(), (n / 2) * ((n + 1) / 2));

    auto bestSolution = solver.getBestSolution();

    int numEdgesInCut = 0;
    for (uint i = 0; i < n; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (uint j = i + 1; j < n; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(g.hasEdge(i, j));
                numEdgesInCut++;
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, (n / 2) * ((n + 1) / 2));
}

TEST(SpecialStructuresSolver, largeApplicableKhalfhalfWeighted) {
    double posEdgeWeight = 1.5;
    uint n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<uint> diss(0, n);

    NetworKit::Graph g(n, true);

    int numMissingEdges = 0;
    for (uint i = 0; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            if (i < n / 2 && j >= n / 2) {
                g.addEdge(i, j, posEdgeWeight);
            } else {
                if (diss(rng) % 2 == 0)
                    g.addEdge(i, j, posEdgeWeight);
                else {
                    if (diss(rng) % 2 == 0)
                        g.addEdge(i, j, -7);
                    numMissingEdges++;
                }
            }
        }
    }

    sms::SpecialStructuresSolver solver(g);
    EXPECT_TRUE(solver.applicable());
    solver.run();

    EXPECT_EQ(solver.bestSolutionValue(), static_cast<double>((n / 2) * ((n + 1) / 2)) * posEdgeWeight);

    auto bestSolution = solver.getBestSolution();

    int numEdgesInCut = 0;
    double cutValue = 0;
    for (uint i = 0; i < n; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (uint j = i + 1; j < n; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(g.hasEdge(i, j));
                numEdgesInCut++;
                cutValue += g.weight(i, j);
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, (n / 2) * ((n + 1) / 2));
    EXPECT_EQ(cutValue, static_cast<double>((n / 2) * ((n + 1) / 2)) * posEdgeWeight);
}

TEST(SpecialStructuresSolver, largeNonApplicableTest) {
    unsigned int n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned int> diss(0, n);

    NetworKit::Graph g(n);

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            if (i <= n / 10 || j != i + 1) {
                if (diss(rng) % 4 != 0)
                    g.addEdge(i, j);
            }
        }
    }

    NetworKit::Graph gCop = g;

    sms::SpecialStructuresSolver solver(g);
    EXPECT_FALSE(solver.applicable());
    EXPECT_ANY_THROW(solver.run());

    EXPECT_EQ(g.numberOfNodes(), gCop.numberOfNodes());
    EXPECT_EQ(g.numberOfEdges(), gCop.numberOfEdges());
}

TEST(SpecialStructuresSolver, supremePositiveEdgeTestSmall) {
    uint n = 5;

    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, -2);
    g.addEdge(0, 3, -17);
    g.addEdge(0, 4, 4);
    g.addEdge(1, 4, -2);
    g.addEdge(2, 4, -1);

    sms::SpecialStructuresSolver solver(g);

    EXPECT_TRUE(solver.applicable());
    solver.run();
    EXPECT_EQ(solver.bestSolutionValue(), 2);

    auto bestSolution = solver.getBestSolution();

    double cutValue = 0;
    for (uint i = 0; i < n; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (uint j = i + 1; j < n; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                if (g.hasEdge(i, j))
                    cutValue += g.weight(i, j);
            }
        }
    }

    EXPECT_EQ(cutValue, solver.bestSolutionValue());
}

TEST(SpecialStructuresSolver, posNegBipartiteTestLarge) {
    uint n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned int> diss(0, n - 1);

    std::vector<uint> nodePermutation(n, 0);
    for (uint i = 0; i < n; i++)
        nodePermutation[i] = i;

    std::shuffle(nodePermutation.begin(), nodePermutation.end(), rng);

    double resultOpt = 0;

    NetworKit::Graph g(n, true);

    for (uint i = 0; i < 10 * n; i++) {
        uint v = diss(rng), w = diss(rng);
        if (v != w && !g.hasEdge(nodePermutation[v], nodePermutation[w])) {
            if ((v < n / 2 && w < n / 2) || (v >= n / 2 && w >= n / 2)) {
                g.addEdge(nodePermutation[v], nodePermutation[w], -1.0 * static_cast<double>(diss(rng)) - 1.0);
            } else {
                double weight = static_cast<double>(diss(rng)) + 1;
                g.addEdge(nodePermutation[v], nodePermutation[w], weight);
                resultOpt += weight;
            }
        }
    }

    sms::SpecialStructuresSolver solver(g);
    EXPECT_TRUE(solver.applicable());
    solver.run();

    EXPECT_EQ(solver.bestSolutionValue(), resultOpt);

    auto bestSolution = solver.getBestSolution();

    double cutValue = 0;
    for (uint i = 0; i < n; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (uint j = i + 1; j < n; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                if (g.hasEdge(i, j))
                    cutValue += g.weight(i, j);
            }
        }
    }

    EXPECT_EQ(cutValue, solver.bestSolutionValue());
}