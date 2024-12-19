#include <gtest/gtest.h>

#include "sms/auxiliary/active_elements.hpp"

TEST(activeElementsStack, empty) {
    auto stack = sms::ActiveElementsStack(5u);
    EXPECT_TRUE(stack.empty());
}

TEST(activeElementsStack, fillAndEmpty) {
    auto stack = sms::ActiveElementsStack(10u);
    stack.activate(1u);
    EXPECT_TRUE(stack.isActive(1u));
    stack.activate(3u);
    EXPECT_TRUE(stack.isActive(3u));
    ASSERT_FALSE(stack.empty());
    stack.deactivate(1u);
    ASSERT_FALSE(stack.empty());
    stack.deactivate(3u);
    EXPECT_TRUE(stack.empty());
}

TEST(activeElementsStack, pop) {
    auto stack = sms::ActiveElementsStack(10u);
    stack.activate(1u);
    EXPECT_TRUE(stack.isActive(1u));
    stack.activate(3u);
    EXPECT_TRUE(stack.isActive(3u));
    ASSERT_FALSE(stack.empty());
    stack.deactivate(1u);
    ASSERT_FALSE(stack.empty());

    unsigned int next = stack.popBack();
    ASSERT_EQ(next, 3u);
}