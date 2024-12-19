#include <gtest/gtest.h>

#include "sms/instance/ising.hpp"

TEST(Ising, createClass1) {
    NetworKit::Graph g(0);
    sms::Ising ising(g);
    ASSERT_EQ(ising.getNumberOfSpins(), 0);
    ASSERT_EQ(ising.getNumberOfInteractions(), 0);
    ASSERT_EQ(ising.getScalingFactor(), 1);
    ASSERT_EQ(ising.getOffset(), 0);
}

TEST(Ising, createClass2) {
    NetworKit::Graph g(3);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);

    sms::Ising ising(g);
    ASSERT_EQ(ising.getNumberOfSpins(), 3);
    ASSERT_EQ(ising.getNumberOfInteractions(), 3);
    ASSERT_EQ(ising.getScalingFactor(), 1);
    ASSERT_EQ(ising.getOffset(), 0);
    ASSERT_EQ(ising.getSolutionValue({1, 1, 1}), 3);
    ASSERT_EQ(ising.getSolutionValue({1, -1, 1}), -1);
    ASSERT_EQ(ising.getSolutionValue({-1, -1, -1}), 3);

    auto info = ising.getInstanceInformation();

    ASSERT_EQ(info["num spins"].get<int>(), 3);
    ASSERT_EQ(info["num interactions"].get<int>(), 3);
    ASSERT_EQ(info["scaling factor"].get<double>(), 1);
    ASSERT_EQ(info["offset"].get<double>(), 0);
}
