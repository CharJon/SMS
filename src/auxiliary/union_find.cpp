#include "sms/auxiliary/union_find.hpp"

#include <cassert>

sms::UnionFind::UnionFind(std::size_t n) : parent_(n), size_(n, 1) {
    for (std::size_t i = 0; i < n; ++i) {
        parent_[i] = i;
    }
}

int sms::UnionFind::numElements() {
    assert(parent_.size() == size_.size());
    return parent_.size();
}

std::size_t sms::UnionFind::find(std::size_t x) {
    assert(x < static_cast<std::size_t>(numElements()));
    // Representative of sets are their own parents
    while (x != parent_[x]) {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return x;
}

void sms::UnionFind::unite(std::size_t x, std::size_t y) {
    auto xRoot = find(x);
    auto yRoot = find(y);

    if (xRoot == yRoot) {
        return;
    }

    if (size_[xRoot] < size_[yRoot]) {
        std::swap(xRoot, yRoot);
    }

    parent_[yRoot] = xRoot;
    size_[xRoot] += size_[yRoot];
}

bool sms::UnionFind::inSameSet(std::size_t x, std::size_t y) {
    return find(x) == find(y);
}
