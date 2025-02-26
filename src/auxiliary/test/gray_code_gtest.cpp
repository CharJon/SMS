#include <random>

#include <gtest/gtest.h>

#include "sms/auxiliary/gray_code.hpp"

TEST(GrayCode, BinaryToGray) {
    EXPECT_EQ(binaryToGray(0u), 0u);
    EXPECT_EQ(binaryToGray(1u), 1u);
    EXPECT_EQ(binaryToGray(2u), 3u);
    EXPECT_EQ(binaryToGray(3u), 2u);
}

TEST(GrayCode, Diff) {
    EXPECT_EQ(diffPosition(0u, 1u), 0);
    EXPECT_EQ(diffPosition(64u, 65u), 0);
    EXPECT_EQ(diffPosition(15u, 7u), 3);
}
