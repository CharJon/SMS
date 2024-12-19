#include <gtest/gtest.h>

#include "networkit/io/EdgeListReader.hpp"

#include "sms/instance/convert.hpp"
#include "sms/instance/frustration_index.hpp"
#include "sms/instance/qubo.hpp"
#include "sms/io/io.hpp"
#include "sms/io/ising_reader.hpp"
#include "sms/solver/enumeration_solver.hpp"

NetworKit::Graph readGraph(const std::string &path) {
    auto er = NetworKit::EdgeListReader(' ', 0, "#");
    auto g = er.read(path);
    assert(g.isWeighted());
    assert(!g.isDirected());
    return g;
}

TEST(MaxCutToQuBO, SymmetryTransformSquare) {
    NetworKit::Graph g = readGraph("test/data/square.wel");
    const int dim = 4;

    sms::MaxCut mc(g);
    auto qubo = maxCutToQUBO(mc);

    double solutionVector[dim] = {1, 0, 0, 1};
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
    ASSERT_EQ(qubo.getOriginalSolutionValue(solutionVector), 1.1);

    solutionVector[0] = 0;
    solutionVector[1] = 1;
    solutionVector[2] = 1;
    solutionVector[3] = 0;
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
}

TEST(MaxCutToQuBO, SymmetryTransformSquareWithOffset) {
    NetworKit::Graph g = readGraph("test/data/square.wel");
    const int dim = 4;

    sms::MaxCut mc(g);
    mc.setOffset(10);
    auto qubo = maxCutToQUBO(mc);

    double solutionVector[dim] = {1, 0, 0, 1};
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
    ASSERT_EQ(qubo.getOriginalSolutionValue(solutionVector), 11.1);

    solutionVector[0] = 0;
    solutionVector[1] = 1;
    solutionVector[2] = 1;
    solutionVector[3] = 0;
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
}

TEST(MaxCutToQuBO, NoSymmetryTransformSquare) {
    NetworKit::Graph g = readGraph("test/data/square.wel");
    const int dim = 3;

    sms::MaxCut mc(g);
    auto qubo = sms::maxCutToQUBOrooted(mc);

    double solutionVector[dim] = {0, 1, 1};
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
    ASSERT_EQ(qubo.getOriginalSolutionValue(solutionVector), 1.1);
}

TEST(MaxCutToQuBO, NoSymmetryTransformSquareWithOffset) {
    NetworKit::Graph g = readGraph("test/data/square.wel");
    const int dim = 3;

    sms::MaxCut mc(g);
    mc.setOffset(10);
    auto qubo = sms::maxCutToQUBOrooted(mc);

    double solutionVector[dim] = {0, 1, 1};
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -1.1);
    ASSERT_EQ(qubo.getOriginalSolutionValue(solutionVector), 11.1);
}

TEST(MaxCutToQuBO, MediumTest) {
    NetworKit::Graph g = readGraph("test/data/tri_square_tri.wel");
    const int dim = 6;

    sms::MaxCut mc(g);
    auto qubo = maxCutToQUBO(mc);

    double solutionVector[dim] = {1, 1, 0, 1, 1, 0};
    ASSERT_EQ(qubo.getSolutionValue(solutionVector), -11.5);
}

TEST(QUBOToMaxCut, BasicTest) {
    auto qubo = sms::bqFileToQubo("test/data/square.bq");

    auto mc = sms::QUBOToMaxCut(qubo);

    ASSERT_EQ(mc.getNumberOfVertices(), 5);
    ASSERT_EQ(mc.getNumberOfEdges(), 8);

    auto solver = sms::EnumerationSolver(mc.getGraph());
    solver.run();
    auto solution = solver.getBestSolution();

    std::vector<uint8_t> sol_0_1(5, 0);

    for (int i = 0; i < 5; i++) {
        if (solution[i]) {
            sol_0_1[i] = 1;
        }
    }

    ASSERT_EQ(mc.getScalingFactor() * mc.getSolutionValue(sol_0_1), -2.0);
}

TEST(QUBOToMaxCut, BasicTest2) {
    auto qubo = sms::bqFileToQubo("test/data/three.bq");

    auto mc = sms::QUBOToMaxCut(qubo);

    ASSERT_EQ(mc.getNumberOfVertices(), 4);
    ASSERT_EQ(mc.getNumberOfEdges(), 5);

    auto solver = sms::EnumerationSolver(mc.getGraph());
    solver.run();
    auto solution = solver.getBestSolution();

    std::vector<uint8_t> sol_0_1(4, 0);

    for (int i = 0; i < 4; i++) {
        if (solution[i]) {
            sol_0_1[i] = 1;
        }
    }

    ASSERT_EQ(mc.getScalingFactor() * mc.getSolutionValue(sol_0_1), -4.0);
}

TEST(FrustrationIndexToMaxCut, SmallTest) {
    NetworKit::Graph g(4, true);

    g.addEdge(0, 1, -1);
    g.addEdge(1, 3, -1);
    g.addEdge(3, 2, -1);
    g.addEdge(2, 0, -1);
    g.addEdge(2, 1, 1);

    sms::FrustrationIndex frust(g);

    auto mc = sms::frustrationIndexToMaxCut(frust);

    std::vector<uint8_t> solutionVector(4, 0);
    solutionVector[1] = 1;
    solutionVector[2] = 1;

    ASSERT_EQ(frust.getSolutionValue(solutionVector), 0);
    ASSERT_EQ(mc.getSolutionValue(solutionVector), 0);

    solutionVector[3] = 1;
    ASSERT_EQ(frust.getSolutionValue(solutionVector), mc.getScalingFactor() * mc.getSolutionValue(solutionVector));

    solutionVector[0] = 1;
    ASSERT_EQ(frust.getSolutionValue(solutionVector), mc.getScalingFactor() * mc.getSolutionValue(solutionVector));

    solutionVector[2] = 0;
    ASSERT_EQ(frust.getSolutionValue(solutionVector), mc.getScalingFactor() * mc.getSolutionValue(solutionVector));

    solutionVector[1] = 0;
    ASSERT_EQ(frust.getSolutionValue(solutionVector), mc.getScalingFactor() * mc.getSolutionValue(solutionVector));
}

TEST(FrustrationIndexToMaxCut, SolvingTest) {
    NetworKit::Graph g(4, true);

    g.addEdge(0, 1, -1);
    g.addEdge(1, 3, -1);
    g.addEdge(3, 2, -1);
    g.addEdge(2, 0, -1);
    g.addEdge(2, 1, 1);

    sms::FrustrationIndex frust(g);

    auto mc = sms::frustrationIndexToMaxCut(frust);
    EXPECT_EQ(mc.getOffset(), -4);
    EXPECT_EQ(mc.getScalingFactor() * mc.getSolutionValue({0, 0, 0, 0, 0}), 4.0);

    auto solver = sms::EnumerationSolver(mc.getGraph());
    solver.run();
    auto solution = solver.getBestSolution();

    std::vector<uint8_t> sol_0_1(4, 0);

    for (int i = 0; i < 4; i++) {
        if (solution[i]) {
            sol_0_1[i] = 1;
        }
    }

    ASSERT_EQ(mc.getScalingFactor() * mc.getSolutionValue(sol_0_1), 0.0);
}

TEST(IsingToMaxCut, SimpleTest) {
    sms::SgParser parser("test/data/simpleGrid.gsg");

    auto ising = parser.getInstance();

    auto mc = sms::isingToMaxCut(ising);
    EXPECT_EQ(mc.getOffset(), 18.0);
    EXPECT_EQ(mc.getSolutionValue({0, 0, 0, 0, 0, 0, 0, 0, 0}), 18.0);

    auto solver = sms::EnumerationSolver(mc.getGraph());
    solver.run();
    auto solution = solver.getBestSolution();

    std::vector<uint8_t> sol_0_1(9, 0);

    for (int i = 0; i < 9; i++) {
        if (solution[i]) {
            sol_0_1[i] = 1;
        }
    }

    ASSERT_EQ(mc.getSolutionValue(sol_0_1), 18.0);
}

TEST(IsingToMaxCut, MediumTest) {
    sms::GsgParser parser("test/data/t3g7_6666.gsg");

    auto torusInstance = parser.getTorus();
    ASSERT_TRUE(torusInstance.has_value());
    if (torusInstance.has_value()) {
        auto ising = torusInstance.value();
        auto mc = sms::isingToMaxCut(ising);

        std::vector<uint8_t> sol(343, 0);

        // solution and value reported by McGroundState
        for (auto i :
             {1,   2,   5,   6,   7,   8,   9,   11,  13,  14,  15,  18,  24,  25,  26,  28,  32,  34,  38,  39,  40,
              41,  42,  43,  44,  46,  47,  49,  50,  51,  54,  58,  61,  64,  72,  74,  75,  77,  81,  85,  88,  92,
              93,  94,  96,  100, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 113, 115, 116, 120, 122, 123, 127,
              128, 130, 131, 137, 138, 148, 149, 150, 152, 153, 154, 160, 161, 163, 169, 170, 171, 172, 174, 177, 180,
              182, 183, 187, 188, 191, 192, 194, 195, 197, 199, 201, 205, 212, 213, 215, 218, 220, 222, 224, 225, 227,
              232, 233, 235, 236, 238, 239, 240, 242, 243, 244, 246, 247, 249, 253, 254, 256, 260, 261, 262, 265, 266,
              272, 274, 275, 276, 279, 280, 281, 282, 283, 285, 288, 289, 293, 294, 295, 302, 303, 304, 305, 307, 309,
              310, 311, 314, 315, 320, 321, 325, 327, 328, 329, 331, 332, 340, 342})
            sol[i - 1] = 1;

        ASSERT_EQ(mc.getSolutionValue(sol), 59074079.0);
    }
}

TEST(graphToTorusIsing, 2DReject) {
    auto g = mcFileToGraph("test/data/nonGrid2d.mc");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_FALSE(converter.isGrid());
}

TEST(graphToTorusIsing, 2DAccept) {
    auto g = mcFileToGraph("test/data/grid2d.mc");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_TRUE(converter.isGrid());

    sms::torusIsing res = converter.getTorusInstance();

    ASSERT_TRUE(res.getGraph().isWeighted());
    ASSERT_TRUE(res.isConsistent());
}

TEST(graphToTorusIsing, 3DAccept) {
    auto g = mcFileToGraph("test/data/t3g7_6666.gsg");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_TRUE(converter.isGrid());

    sms::torusIsing res = converter.getTorusInstance();

    ASSERT_TRUE(res.getGraph().isWeighted());
    ASSERT_TRUE(res.isConsistent());
}

TEST(graphToTorusIsing, 3DAcceptRelabeled) {
    auto g = mcFileToGraph("test/data/t3g7_6666.mc");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_TRUE(converter.isGrid());

    sms::torusIsing res = converter.getTorusInstance();

    ASSERT_TRUE(res.getGraph().isWeighted());
    ASSERT_TRUE(res.isConsistent());
}

TEST(graphToTorusIsing, 3DAcceptRandomRelabeled) {
    auto g = mcFileToGraph("test/data/t3g7_6666_random_label.mc");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_TRUE(converter.isGrid());

    sms::torusIsing res = converter.getTorusInstance();

    ASSERT_TRUE(res.getGraph().isWeighted());
    ASSERT_TRUE(res.isConsistent());
}

TEST(graphToTorusIsing, 3DReject) {
    auto g = mcFileToGraph("test/data/t3g7_6666_no_grid.mc");
    auto cg = sms::compactGraph(g);
    g = cg.compactGraph;

    sms::graphToTorusIsing converter(g);

    converter.run();

    ASSERT_TRUE(converter.hasRun());
    ASSERT_FALSE(converter.isGrid());
}