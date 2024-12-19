#ifndef SMS_UNION_FIND_HPP
#define SMS_UNION_FIND_HPP

#include "vector"

namespace sms {

class UnionFind {

public:
    explicit UnionFind(std::size_t n);

    int numElements();

    std::size_t find(std::size_t x);

    void unite(std::size_t x, std::size_t y);

    bool inSameSet(std::size_t x, std::size_t y);

private:
    std::vector<std::size_t> parent_;
    std::vector<int> size_;
};

} // namespace sms

#endif // SMS_UNION_FIND_HPP
