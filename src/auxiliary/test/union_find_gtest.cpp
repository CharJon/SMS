#include <gtest/gtest.h>

#include "sms/auxiliary/union_find.hpp"

TEST(unionFind, init) {
    int numElements = 5;
    auto uf = sms::UnionFind(numElements);
    ASSERT_EQ(uf.numElements(), numElements);

    for (int i = 0; i < numElements; ++i) {
        EXPECT_EQ(uf.find(i), i);
    }
}

TEST(unionFind, fillAndEmpty) {
    auto uf = sms::UnionFind(10);

    EXPECT_EQ(uf.find(1), 1);
    uf.unite(1, 2);
    uf.unite(1, 3);
    EXPECT_EQ(uf.find(2), uf.find(3));
    uf.unite(4, 1);
    uf.unite(5, 6);
    uf.unite(5, 7);
    uf.unite(7, 1);
    EXPECT_TRUE(uf.inSameSet(1, 7));
}
