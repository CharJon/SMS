#include <gtest/gtest.h>

#include "sms/io/io.hpp"
#include "sms/io/mc_read_write.hpp"

TEST(ReadWriteMC, test) {
    NetworKit::Graph g(5, true);
    g.addEdge(0, 0);
    auto tmpFile = TmpFile();
    ASSERT_ANY_THROW(sms::writeToMcFile(g, tmpFile.getFilename()));
}

TEST(ReadWriteMC, SimpleSaveTest) {
    // read
    auto orig = mcFileToGraph("test/data/square.mc");
    // write
    auto tmpFile = TmpFile();
    sms::writeToMcFile(orig, tmpFile.getFilename());
    // read again
    auto g = mcFileToGraph(tmpFile.getFilename());
    // graph is not compacted, number of vertices increases by 1
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 4);
    for (auto e : orig.edgeWeightRange()) {
        ASSERT_TRUE(g.hasEdge(e.u + 1, e.v + 1));
        ASSERT_EQ(g.weight(e.u + 1, e.v + 1), e.weight);
    }
}

TEST(ReadWriteMC, ComplexSaveTest) {
    auto orig = mcFileToGraph("test/data/mannino_k487b.mc");
    auto tmpFile = TmpFile();
    sms::writeToMcFile(orig, tmpFile.getFilename());

    auto g = mcFileToGraph(tmpFile.getFilename());
    ASSERT_EQ(g.numberOfNodes(), orig.numberOfNodes() + 1);
    ASSERT_EQ(g.numberOfEdges(), orig.numberOfEdges());

    for (auto e : orig.edgeWeightRange()) {
        ASSERT_TRUE(g.hasEdge(e.u + 1, e.v + 1));
        ASSERT_EQ(g.weight(e.u + 1, e.v + 1), e.weight);
    }
}

TEST(ReadWriteMC, ErrorSaveTestSelfloop) {
    NetworKit::Graph orig(2, true);
    orig.addEdge(1, 1, 1);
    auto tmpFile = TmpFile();
    ASSERT_ANY_THROW(sms::writeToMcFile(orig, tmpFile.getFilename()));
}

TEST(ReadWriteMC, DoubleEdgeTest) {
    NetworKit::Graph orig(3, true);
    orig.addEdge(1, 2, 1);
    orig.addEdge(1, 2, 1);
    auto tmpFile = TmpFile();
    ASSERT_ANY_THROW(sms::writeToMcFile(orig, tmpFile.getFilename()));
}

TEST(ReadWriteMC, DoubleEdgeTest2) {
    NetworKit::Graph orig(3, true);
    orig.addEdge(1, 2, 1);
    orig.addEdge(2, 1, 1);
    auto tmpFile = TmpFile();
    ASSERT_ANY_THROW(sms::writeToMcFile(orig, tmpFile.getFilename()));
}

TEST(ReadWriteMC, DoubleEdgeTest3) {
    NetworKit::Graph orig(3, true);
    orig.addEdge(1, 2, 1);
    orig.addEdge(2, 1, 0);
    auto tmpFile = TmpFile();
    ASSERT_ANY_THROW(sms::writeToMcFile(orig, tmpFile.getFilename()));
}
