#include <gtest/gtest.h>

#include "sms/conshdlr/cliques.hpp"

TEST(CliqueSeparator, fiveViolated) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, 0.4);
    g.addEdge(0, 2, 0.4);
    g.addEdge(0, 3, 0.4);
    g.addEdge(0, 4, 0.4);
    g.addEdge(1, 2, 0.4);
    g.addEdge(1, 3, 0.4);
    g.addEdge(1, 4, 0.4);
    g.addEdge(2, 3, 0.4);
    g.addEdge(2, 4, 0.4);
    g.addEdge(3, 4, 0.4);

    auto cs = sms::CliqueSeparator(g);
    cs.addClique({0, 1, 2, 3, 4});

    cs.updateWeights(0, 1, 2. / 3);
    cs.updateWeights(0, 2, 2. / 3);
    cs.updateWeights(0, 3, 2. / 3);
    cs.updateWeights(0, 4, 2. / 3);
    cs.updateWeights(1, 2, 2. / 3);
    cs.updateWeights(1, 3, 2. / 3);
    cs.updateWeights(1, 4, 2. / 3);
    cs.updateWeights(2, 3, 2. / 3);
    cs.updateWeights(2, 4, 2. / 3);
    cs.updateWeights(3, 4, 2. / 3);
    auto checkViolation = cs.checkViolation(cs.getClique(0));
    auto res = checkViolation[0];
    EXPECT_EQ(res.first.size() + res.second.size(), 5);
    EXPECT_TRUE(!res.first.empty() || !res.second.empty());
}

TEST(CliqueSeparator, fiveNoViolated) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 1, 0.4);
    g.addEdge(0, 2, 0.4);
    g.addEdge(0, 3, 0.4);
    g.addEdge(0, 4, 0.4);
    g.addEdge(1, 2, 0.4);
    g.addEdge(1, 3, 0.4);
    g.addEdge(1, 4, 0.4);
    g.addEdge(2, 3, 0.4);
    g.addEdge(2, 4, 0.4);
    g.addEdge(3, 4, 0.4);

    auto cs = sms::CliqueSeparator(g);
    cs.addClique({0, 1, 2, 3, 4});

    cs.updateWeights(0, 1, 0);
    cs.updateWeights(0, 2, 0);
    cs.updateWeights(0, 3, 0);
    cs.updateWeights(0, 4, 0);
    cs.updateWeights(1, 2, 0);
    cs.updateWeights(1, 3, 0);
    cs.updateWeights(1, 4, 0);
    cs.updateWeights(2, 3, 0);
    cs.updateWeights(2, 4, 0);
    cs.updateWeights(3, 4, 0);
    auto checkViolation = cs.checkViolation(cs.getClique(0));
    EXPECT_TRUE(checkViolation.empty());
}

TEST(CliqueSeparator, bigOddClique) {
    NetworKit::Graph g(19, true);
    std::vector<NetworKit::node> clique = {};
    for (unsigned int i = 0; i < 19; i++) {
        clique.push_back(i);
        for (unsigned int j = i + 1; j < 19; j++) {
            g.addEdge(i, j, 0.1);
        }
    }

    auto cs = sms::CliqueSeparator(g);
    cs.addClique(clique);

    for (int i = 0; i < 15; ++i) {
        for (int j = i + 1; j < 15; ++j) {
            cs.updateWeights(i, j, 2. / 3);
        }
    }
    auto checkViolation = cs.checkViolation(cs.getClique(0));
    auto res = checkViolation[0];
    // EXPECT_EQ(res.first.size() + res.second.size(), 5);
    EXPECT_TRUE(!res.first.empty() || !res.second.empty());
}

TEST(CliqueSeparator, bigClique1) {
    NetworKit::Graph g(25, true);
    std::vector<NetworKit::node> clique = {};
    for (unsigned int i = 0; i < 25; i++) {
        clique.push_back(i);
        for (unsigned int j = i + 1; j < 25; j++) {
            g.addEdge(i, j, 0.5);
        }
    }

    auto cs = sms::CliqueSeparator(g);
    cs.addClique(clique);
    cs.updateWeights(1, 2, 1);
    cs.updateWeights(1, 3, 2. / 3);
    cs.updateWeights(1, 4, 2. / 3);
    cs.updateWeights(2, 4, 2. / 3);
    cs.updateWeights(3, 4, 2. / 3);
    auto checkViolation = cs.checkViolation(cs.getClique(0));
    auto res = checkViolation[0];
    EXPECT_EQ(res.first.size() + res.second.size(), 5);
    EXPECT_TRUE(!res.first.empty() || !res.second.empty());
}

TEST(CliqueSeparator, bigClique2) {
    NetworKit::Graph g(25, true);
    std::vector<NetworKit::node> clique = {};
    for (unsigned int i = 0; i < 25; i++) {
        clique.push_back(i);
        for (unsigned int j = i + 1; j < 25; j++) {
            g.addEdge(i, j, 0.5);
        }
    }

    auto cs = sms::CliqueSeparator(g);
    cs.addClique(clique);
    cs.updateWeights(1, 2, 1);
    cs.updateWeights(1, 3, 0.1);
    cs.updateWeights(0, 2, 0.8);
    cs.updateWeights(0, 3, 0.4);
    cs.updateWeights(10, 9, 1);
    cs.updateWeights(11, 9, 0.4);
    auto checkViolation = cs.checkViolation(cs.getClique(0));
    auto res = checkViolation[0];
    EXPECT_TRUE(!res.first.empty() || !res.second.empty());
}