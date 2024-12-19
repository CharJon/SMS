#include <gtest/gtest.h>

#include "sms/auxiliary/chrono_io.hpp"

TEST(chronoIO, BasicTest) {
    auto s = std::chrono::seconds(1);

    std::stringstream ss;
    ss << s;
    ASSERT_EQ(ss.str(), "1s");
}
