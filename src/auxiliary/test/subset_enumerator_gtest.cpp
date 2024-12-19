#include <gtest/gtest.h>

#include "sms/auxiliary/subset_enumerator.hpp"

int cntOneBits(uint64_t n) {
    int cnt = 0;
    while (n) {
        cnt += n & 1;
        n >>= 1;
    }
    return cnt;
}

TEST(subsetEnumerator, fiveThree) {
    auto se = sms::SubsetEnumerator(5, 3);
    int cnt = 0;
    while (!se.done()) {
        cnt++;
        auto next = se.next();
        EXPECT_EQ(cntOneBits(next), 3);
    }
    ASSERT_EQ(cnt, 10);
}

TEST(subsetEnumerator, elevenFive) {
    auto se = sms::SubsetEnumerator(11, 5);
    int cnt = 0;
    while (!se.done()) {
        cnt++;
        auto next = se.next();
        EXPECT_EQ(cntOneBits(next), 5);
    }
    ASSERT_EQ(cnt, 462);
}