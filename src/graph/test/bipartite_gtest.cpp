#include <gtest/gtest.h>

#include "sms/graph/bipartite.hpp"

TEST(Bipartite, SimpleGraph1) {
    auto g = NetworKit::Graph(5);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(2, 4);
    g.removeNode(3);

    sms::Bipartite bipar(g);
    bipar.run();
    ASSERT_TRUE(bipar.isBipartiteGraph());
    const auto &partition = bipar.getPartition();
    EXPECT_NE(partition[0], partition[1]);
    EXPECT_NE(partition[0], partition[2]);
    EXPECT_EQ(partition[0], partition[4]);
}

TEST(Bipartite, SimpleGraphNegativeWeights) {
    auto g = NetworKit::Graph(3, true, false);

    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(1, 2, -1);

    sms::Bipartite bipar(g);
    bipar.run();
    EXPECT_FALSE(bipar.isBipartiteGraph());
    bipar.run(true);
    EXPECT_TRUE(bipar.isBipartiteGraph());
}
