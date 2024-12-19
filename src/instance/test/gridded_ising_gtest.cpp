#include <gtest/gtest.h>

#include "sms/instance/gridded_ising.hpp"

TEST(griddedIsing, createClass1) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);
    ASSERT_EQ(gIsing.getNumberOfSpins(), 0);
    ASSERT_EQ(gIsing.getNumberOfInteractions(), 0);
    ASSERT_EQ(gIsing.getScalingFactor(), 1);
    ASSERT_EQ(gIsing.getOffset(), 0);
    ASSERT_FALSE(gIsing.isConsistent());

    auto info = gIsing.getInstanceInformation();

    ASSERT_EQ(info["num spins"].get<int>(), 0);
    ASSERT_EQ(info["num interactions"].get<int>(), 0);
    ASSERT_EQ(info["scaling factor"].get<double>(), 1);
    ASSERT_EQ(info["offset"].get<double>(), 0);
    ASSERT_EQ(info["dimension"].get<int>(), 2);
    ASSERT_EQ(info["grid length"].get<int>(), 3);
}

TEST(griddedIsing, coordinate2D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);

    ASSERT_EQ(gIsing.getCoordinates(0), sms::gridCoordinates(0, 0, 0));
    ASSERT_EQ(gIsing.getCoordinates(3), sms::gridCoordinates(0, 1, 0));
    ASSERT_EQ(gIsing.getCoordinates(7), sms::gridCoordinates(0, 2, 1));
    ASSERT_EQ(gIsing.getCoordinates(8), sms::gridCoordinates(0, 2, 2));

    ASSERT_EQ(gIsing.getNodeAt({0, 0, 0}), 0);
    ASSERT_EQ(gIsing.getNodeAt({0, 1, 0}), 3);
    ASSERT_EQ(gIsing.getNodeAt({0, 2, 1}), 7);
    ASSERT_EQ(gIsing.getNodeAt({0, 2, 2}), 8);
}

TEST(griddedIsing, coordinate3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);

    ASSERT_EQ(gIsing.getCoordinates(0), sms::gridCoordinates(0, 0, 0));
    ASSERT_EQ(gIsing.getCoordinates(3), sms::gridCoordinates(0, 1, 0));
    ASSERT_EQ(gIsing.getCoordinates(7), sms::gridCoordinates(0, 2, 1));
    ASSERT_EQ(gIsing.getCoordinates(8), sms::gridCoordinates(0, 2, 2));

    ASSERT_EQ(gIsing.getCoordinates(9), sms::gridCoordinates(1, 0, 0));
    ASSERT_EQ(gIsing.getCoordinates(12), sms::gridCoordinates(1, 1, 0));
    ASSERT_EQ(gIsing.getCoordinates(16), sms::gridCoordinates(1, 2, 1));
    ASSERT_EQ(gIsing.getCoordinates(17), sms::gridCoordinates(1, 2, 2));

    ASSERT_EQ(gIsing.getCoordinates(18), sms::gridCoordinates(2, 0, 0));
    ASSERT_EQ(gIsing.getCoordinates(21), sms::gridCoordinates(2, 1, 0));
    ASSERT_EQ(gIsing.getCoordinates(25), sms::gridCoordinates(2, 2, 1));
    ASSERT_EQ(gIsing.getCoordinates(26), sms::gridCoordinates(2, 2, 2));

    ASSERT_EQ(gIsing.getNodeAt({0, 0, 0}), 0);
    ASSERT_EQ(gIsing.getNodeAt({0, 1, 0}), 3);
    ASSERT_EQ(gIsing.getNodeAt({0, 2, 1}), 7);
    ASSERT_EQ(gIsing.getNodeAt({0, 2, 2}), 8);

    ASSERT_EQ(gIsing.getNodeAt({1, 0, 0}), 9);
    ASSERT_EQ(gIsing.getNodeAt({1, 1, 0}), 12);
    ASSERT_EQ(gIsing.getNodeAt({1, 2, 1}), 16);
    ASSERT_EQ(gIsing.getNodeAt({1, 2, 2}), 17);

    ASSERT_EQ(gIsing.getNodeAt({2, 0, 0}), 18);
    ASSERT_EQ(gIsing.getNodeAt({2, 1, 0}), 21);
    ASSERT_EQ(gIsing.getNodeAt({2, 2, 1}), 25);
    ASSERT_EQ(gIsing.getNodeAt({2, 2, 2}), 26);
}

TEST(griddedIsing, right2D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);
    ASSERT_EQ(gIsing.right(0), 1);
    ASSERT_EQ(gIsing.right(2), 0);
    ASSERT_EQ(gIsing.right(8), 6);
}

TEST(griddedIsng, right3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.right(0), 1);
    ASSERT_EQ(gIsing.right(2), 0);
    ASSERT_EQ(gIsing.right(8), 6);

    ASSERT_EQ(gIsing.right(9), 10);
    ASSERT_EQ(gIsing.right(11), 9);
    ASSERT_EQ(gIsing.right(17), 15);

    ASSERT_EQ(gIsing.right(18), 19);
    ASSERT_EQ(gIsing.right(20), 18);
    ASSERT_EQ(gIsing.right(26), 24);
}

TEST(griddedIsing, left2D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);
    ASSERT_EQ(gIsing.left(0), 2);
    ASSERT_EQ(gIsing.left(2), 1);
    ASSERT_EQ(gIsing.left(8), 7);
}

TEST(griddedIsng, left3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.left(0), 2);
    ASSERT_EQ(gIsing.left(2), 1);
    ASSERT_EQ(gIsing.left(8), 7);

    ASSERT_EQ(gIsing.left(9), 11);
    ASSERT_EQ(gIsing.left(11), 10);
    ASSERT_EQ(gIsing.left(17), 16);

    ASSERT_EQ(gIsing.left(18), 20);
    ASSERT_EQ(gIsing.left(20), 19);
    ASSERT_EQ(gIsing.left(26), 25);
}

TEST(griddedIsing, up2D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);
    ASSERT_EQ(gIsing.up(0), 6);
    ASSERT_EQ(gIsing.up(2), 8);
    ASSERT_EQ(gIsing.up(8), 5);
}

TEST(griddedIsng, up3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.up(0), 6);
    ASSERT_EQ(gIsing.up(2), 8);
    ASSERT_EQ(gIsing.up(8), 5);

    ASSERT_EQ(gIsing.up(9), 15);
    ASSERT_EQ(gIsing.up(11), 17);
    ASSERT_EQ(gIsing.up(17), 14);

    ASSERT_EQ(gIsing.up(18), 24);
    ASSERT_EQ(gIsing.up(20), 26);
    ASSERT_EQ(gIsing.up(26), 23);
}

TEST(griddedIsing, down2D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 2, 3);
    ASSERT_EQ(gIsing.down(0), 3);
    ASSERT_EQ(gIsing.down(2), 5);
    ASSERT_EQ(gIsing.down(8), 2);
}

TEST(griddedIsng, down3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.down(0), 3);
    ASSERT_EQ(gIsing.down(2), 5);
    ASSERT_EQ(gIsing.down(8), 2);

    ASSERT_EQ(gIsing.down(9), 12);
    ASSERT_EQ(gIsing.down(11), 14);
    ASSERT_EQ(gIsing.down(17), 11);

    ASSERT_EQ(gIsing.down(18), 21);
    ASSERT_EQ(gIsing.down(20), 23);
    ASSERT_EQ(gIsing.down(26), 20);
}

TEST(griddedIsng, back3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.back(0), 9);
    ASSERT_EQ(gIsing.back(2), 11);
    ASSERT_EQ(gIsing.back(8), 17);

    ASSERT_EQ(gIsing.back(9), 18);
    ASSERT_EQ(gIsing.back(11), 20);
    ASSERT_EQ(gIsing.back(17), 26);

    ASSERT_EQ(gIsing.back(18), 0);
    ASSERT_EQ(gIsing.back(20), 2);
    ASSERT_EQ(gIsing.back(26), 8);
}

TEST(griddedIsng, front3D) {
    NetworKit::Graph g(0);
    sms::torusIsing gIsing(g, 3, 3);
    ASSERT_EQ(gIsing.front(0), 18);
    ASSERT_EQ(gIsing.front(2), 20);
    ASSERT_EQ(gIsing.front(8), 26);

    ASSERT_EQ(gIsing.front(9), 0);
    ASSERT_EQ(gIsing.front(11), 2);
    ASSERT_EQ(gIsing.front(17), 8);

    ASSERT_EQ(gIsing.front(18), 9);
    ASSERT_EQ(gIsing.front(20), 11);
    ASSERT_EQ(gIsing.front(26), 17);
}

TEST(computeTorusParams, 2DTori) {
    auto res = sms::computeTorusParams(9, 18);
    ASSERT_TRUE(res.has_value());
    auto [d, l] = res.value();
    ASSERT_EQ(d, 2);
    ASSERT_EQ(l, 3);

    res = sms::computeTorusParams(10, 24);
    ASSERT_FALSE(res.has_value());
}

TEST(computeTorusParams, 3DTori) {
    auto res = sms::computeTorusParams(343, 1029);
    ASSERT_TRUE(res.has_value());

    auto [d, l] = res.value();
    ASSERT_EQ(d, 3);
    ASSERT_EQ(l, 7);

    res = sms::computeTorusParams(1000, 2000);
    ASSERT_FALSE(res.has_value());
}

TEST(computeTorusParams, TorusAmbiguous) {
    auto res = sms::computeTorusParams(36, 72);
    ASSERT_TRUE(res.has_value());
    auto [d, l] = res.value();
    ASSERT_EQ(d, 2);
    ASSERT_EQ(l, 6);

    auto res2 = sms::computeTorusParams(64, 192);
    ASSERT_TRUE(res2.has_value());
    auto [d2, l2] = res2.value();
    ASSERT_EQ(d2, 3);
    ASSERT_EQ(l2, 4);
}

TEST(computeGridParams, 2DGrids) {
    auto res = sms::computeGridParams(9, 12);
    ASSERT_TRUE(res.has_value());
    auto [d, l] = res.value();
    ASSERT_EQ(d, 2);
    ASSERT_EQ(l, 3);

    auto res2 = sms::computeGridParams(1000000, 1998000);
    ASSERT_TRUE(res2.has_value());
    auto [d2, l2] = res2.value();
    ASSERT_EQ(d2, 2);
    ASSERT_EQ(l2, 1000);

    res = sms::computeGridParams(9, 18);
    ASSERT_FALSE(res.has_value());
}

TEST(computeGridParams, 3DGrids) {
    auto res = sms::computeGridParams(27, 54);
    ASSERT_TRUE(res.has_value());
    auto [d, l] = res.value();
    ASSERT_EQ(d, 3);
    ASSERT_EQ(l, 3);

    auto res2 = sms::computeGridParams(1000000, 3000000 - 30000);
    ASSERT_TRUE(res2.has_value());
    auto [d2, l2] = res2.value();
    ASSERT_EQ(d2, 3);
    ASSERT_EQ(l2, 100);
}

TEST(computeGridParams, GridAmbiguous) {
    auto res = sms::computeGridParams(64, 128 - 2 * 8);
    ASSERT_TRUE(res.has_value());
    auto [d, l] = res.value();
    ASSERT_EQ(d, 2);
    ASSERT_EQ(l, 8);

    auto res2 = sms::computeGridParams(64, 192 - 3 * 16);
    ASSERT_TRUE(res2.has_value());
    auto [d2, l2] = res2.value();
    ASSERT_EQ(d2, 3);
    ASSERT_EQ(l2, 4);
}
