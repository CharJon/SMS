#include <gtest/gtest.h>

#include "sms/io/qplib_parser.hpp"

TEST(QPLibToQUBO, ParsingTest) {
    sms::QPLibParser p("test/data/QPLIB_3506.qplib");

    auto i = p.getInstance();

    ASSERT_TRUE(i.isValid());
    EXPECT_EQ(i.getScalingFactor(), -1); // original problem is maximization problem
    EXPECT_EQ(i.getValue(274, 11), 2.0);
    EXPECT_EQ(i.getValue(1, 0), -2.0);
    EXPECT_EQ(i.getValue(426, 12), 2.0);
    EXPECT_EQ(i.getValue(12, 426), 0.0);
    EXPECT_EQ(i.getValue(170, 170), -2.0);
    EXPECT_EQ(i.getValue(251, 251), 2.0);
}
