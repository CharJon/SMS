#include <gtest/gtest.h>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/odd_closed_walk.hpp"

TEST(OddClosedWalk, Init) {
    sms::OddClosedWalk ow(0);
    ow.addEdge(1, false);
    ow.addEdge(2, true);
    ow.addEdge(3, false);
    ow.addEdge(0, false);
}

TEST(OddClosedWalk, isValid) {
    sms::OddClosedWalk ow(0);
    ow.addEdge(1, false);
    ow.addEdge(2, true);
    ow.addEdge(3, true);
    ow.addEdge(0, true);
    ow.addEdge(4, true);
    ow.addEdge(5, true);
    ow.addEdge(6, true);
    ow.addEdge(2, true);
    ow.addEdge(1, true);
    ow.addEdge(2, true);

    ASSERT_FALSE(ow.isValid());

    ow.addEdge(0, false);

    ASSERT_FALSE(ow.isValid());

    ow.addEdge(1, false);
    ow.addEdge(0, true);

    ASSERT_TRUE(ow.isValid());
}

TEST(OddClosedWalk, IsValidGraph) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(1, 0);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 0);
    g.addEdge(0, 3);

    sms::OddClosedWalk ow(0);
    ow.addEdge(1, true);
    ow.addEdge(2, false);
    ow.addEdge(0, true);

    ASSERT_TRUE(ow.isValid());
    ASSERT_FALSE(ow.isValid(g));

    sms::OddClosedWalk ow2(0);

    ow2.addEdge(1, true);
    ow2.addEdge(2, true);
    ow2.addEdge(3, true);
    ow2.addEdge(0, false);
    ow2.addEdge(4, true);
    ow2.addEdge(3, true);
    ow2.addEdge(2, true);
    ow2.addEdge(1, true);
    ow2.addEdge(0, true);

    ASSERT_TRUE(ow2.isValid(g));
}

TEST(OddClosedWalk, Extract) {
    sms::OddClosedWalk ow(0);

    ow.addEdge(1, true);
    ow.addEdge(2, true);
    ow.addEdge(3, true);
    ow.addEdge(0, false);
    ow.addEdge(4, true);
    ow.addEdge(3, true);
    ow.addEdge(2, true);
    ow.addEdge(1, true);
    ow.addEdge(0, true);

    auto ex = ow.extract(3, 6);
    ASSERT_TRUE(ex.isValid());

    ASSERT_DEBUG_DEATH(ow.extract(2, 6), "Assertion");
    ASSERT_DEBUG_DEATH(ow.extract(6, 3), "Assertion");
}

// TEST(OddClosedWalk, splitOnChordInvalidChord) {
//     OddClosedWalk ow(0);
//
//     ow.addEdge(1);
//     ow.addEdge(2);
//     ow.addEdge(3);
//     ow.addEdge(0);
//     ow.addEdge(4);
//     ow.addEdge(3);
//     ow.addEdge(2);
//     ow.addEdge(1);
//     ow.addEdge(0);
//
//     ASSERT_DEBUG_DEATH(ow.splitOnChord(0, 3, true, false), "Assertion");
// }

TEST(OddClosedWalk, splitOnChordSimple) {
    sms::OddClosedWalk ow(0);

    ow.addEdge(1, false);
    ow.addEdge(2, true);
    ow.addEdge(3, true);
    ow.addEdge(4, true);
    ow.addEdge(0, true);

    auto split = ow.splitOnChord(0, 3, false, true);
    ASSERT_EQ(split.getIthNode(0), 0);
    ASSERT_EQ(split.getIthNode(1), 3);
    ASSERT_EQ(split.getIthNode(2), 4);
    ASSERT_EQ(split.getIthNode(3), 0);

    ASSERT_FALSE(split.ithEdgeIsStay(0));
    ASSERT_TRUE(split.isValid());

    auto splitInner = ow.splitOnChord(0, 3, true, false);
    ASSERT_EQ(splitInner.getIthNode(0), 0);
    ASSERT_EQ(splitInner.getIthNode(1), 1);
    ASSERT_EQ(splitInner.getIthNode(2), 2);
    ASSERT_EQ(splitInner.getIthNode(3), 3);
    ASSERT_EQ(splitInner.getIthNode(4), 0);

    ASSERT_TRUE(splitInner.isValid());
}

// TEST(OddClosedWalk, splitOnChordNoCrossEdges) {
//     OddClosedWalk ow(0);
//
//     ow.addEdge(1);
//     ow.addEdge(2);
//     ow.addEdge(3);
//     ow.addEdge(4);
//     ow.addEdge(0);
//
//     ASSERT_DEBUG_DEATH(ow.splitOnChord(0, 3, false, false), "Assertion");
// }

TEST(OddClosedWalk, violatedOddSelectionClosedWalkEvenCross) {
    auto ocw = sms::violatedOddSelectionClosedWalk({1, 2, 3, 4}, {0.95, 0.95, 0.95, 0.6});

    ASSERT_TRUE(ocw.has_value());
    ASSERT_TRUE(ocw->isValid());
    ASSERT_EQ(ocw->size(), 4);
}

TEST(OddClosedWalk, violatedOddSelectionClosedWalkEvenCrossImpossible) {
    auto ocw = sms::violatedOddSelectionClosedWalk({1, 2, 3, 4}, {0.95, 0.95, 0.95, 0.95});

    ASSERT_FALSE(ocw.has_value());
}

TEST(OddClosedWalk, violatedOddSelectionClosedWalkStayToCross) {
    auto ocw = sms::violatedOddSelectionClosedWalk({1, 2, 3, 4}, {0.95, 0.95, 0.1, 0.4});

    ASSERT_TRUE(ocw.has_value());
    ASSERT_TRUE(ocw->isValid());
    ASSERT_EQ(ocw->size(), 4);
}

TEST(OddClosedWalk, violatedOddSelectionClosedWalkAllStay) {
    auto ocw = sms::violatedOddSelectionClosedWalk({1, 2, 3, 4}, {0.1, 0.1, 0.4, 0.1});

    ASSERT_TRUE(ocw.has_value());
    ASSERT_TRUE(ocw->isValid());
    ASSERT_EQ(ocw->size(), 4);
}