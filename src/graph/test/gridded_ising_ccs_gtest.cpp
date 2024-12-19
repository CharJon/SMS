#include <vector>

#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/gridded_ising_ccs.hpp"

TEST(GriddedIsingCCs, 2DTest) {
    NetworKit::Graph g(9);
    sms::torusIsing gIsing(g, 2, 3);

    sms::GriddedIsingChordlessCycles ccs(gIsing);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());

    ASSERT_EQ(ccs.fourHoles.size(), 9);
    ASSERT_EQ(ccs.sixCycles.size(), 0);
    ASSERT_EQ(ccs.eightCycles.size(), 9);

    ASSERT_EQ(ccs.fourHoles[0][0], 0);
    ASSERT_EQ(ccs.fourHoles[0][1], 1);
    ASSERT_EQ(ccs.fourHoles[0][2], 4);
    ASSERT_EQ(ccs.fourHoles[0][3], 3);

    ASSERT_EQ(ccs.eightCycles[0][0], 0);
    ASSERT_EQ(ccs.eightCycles[0][1], 1);
    ASSERT_EQ(ccs.eightCycles[0][2], 2);
    ASSERT_EQ(ccs.eightCycles[0][3], 5);
    ASSERT_EQ(ccs.eightCycles[0][4], 8);
    ASSERT_EQ(ccs.eightCycles[0][5], 7);
    ASSERT_EQ(ccs.eightCycles[0][6], 6);
    ASSERT_EQ(ccs.eightCycles[0][7], 3);
}

TEST(GriddedIsingCCs, 3DTest) {
    NetworKit::Graph g(27);
    sms::torusIsing gIsing(g, 3, 3);

    sms::GriddedIsingChordlessCycles ccs(gIsing);

    ccs.run();

    ASSERT_TRUE(ccs.hasRun());

    ASSERT_EQ(ccs.fourHoles.size(), 3 * 27);
    ASSERT_EQ(ccs.sixCycles.size(), 4 * 27);
    ASSERT_EQ(ccs.eightCycles.size(), 3 * 27);

    ASSERT_EQ(ccs.fourHoles[1][0], 0);
    ASSERT_EQ(ccs.fourHoles[1][1], 9);
    ASSERT_EQ(ccs.fourHoles[1][2], 12);
    ASSERT_EQ(ccs.fourHoles[1][3], 3);

    ASSERT_EQ(ccs.eightCycles[1][0], 0);
    ASSERT_EQ(ccs.eightCycles[1][1], 9);
    ASSERT_EQ(ccs.eightCycles[1][2], 18);
    ASSERT_EQ(ccs.eightCycles[1][3], 21);
    ASSERT_EQ(ccs.eightCycles[1][4], 24);
    ASSERT_EQ(ccs.eightCycles[1][5], 15);
    ASSERT_EQ(ccs.eightCycles[1][6], 6);
    ASSERT_EQ(ccs.eightCycles[1][7], 3);

    ASSERT_EQ(ccs.sixCycles[0][0], 0);
    ASSERT_EQ(ccs.sixCycles[0][1], 9);
    ASSERT_EQ(ccs.sixCycles[0][2], 10);
    ASSERT_EQ(ccs.sixCycles[0][3], 13);
    ASSERT_EQ(ccs.sixCycles[0][4], 4);
    ASSERT_EQ(ccs.sixCycles[0][5], 3);

    ASSERT_EQ(ccs.sixCycles[1][0], 0);
    ASSERT_EQ(ccs.sixCycles[1][1], 1);
    ASSERT_EQ(ccs.sixCycles[1][2], 4);
    ASSERT_EQ(ccs.sixCycles[1][3], 13);
    ASSERT_EQ(ccs.sixCycles[1][4], 12);
    ASSERT_EQ(ccs.sixCycles[1][5], 9);

    ASSERT_EQ(ccs.sixCycles[2][0], 0);
    ASSERT_EQ(ccs.sixCycles[2][1], 3);
    ASSERT_EQ(ccs.sixCycles[2][2], 12);
    ASSERT_EQ(ccs.sixCycles[2][3], 13);
    ASSERT_EQ(ccs.sixCycles[2][4], 10);
    ASSERT_EQ(ccs.sixCycles[2][5], 1);

    ASSERT_EQ(ccs.sixCycles[3][0], 0);
    ASSERT_EQ(ccs.sixCycles[3][1], 9);
    ASSERT_EQ(ccs.sixCycles[3][2], 11);
    ASSERT_EQ(ccs.sixCycles[3][3], 14);
    ASSERT_EQ(ccs.sixCycles[3][4], 5);
    ASSERT_EQ(ccs.sixCycles[3][5], 3);
}
