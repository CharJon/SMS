#ifndef SMS_GRAY_CODE_HPP
#define SMS_GRAY_CODE_HPP

#include <bit>
#include <cassert>
#include <concepts>

template <typename T>
requires std::unsigned_integral<T>
T binaryToGray(T num) {
    return num ^ (num >> 1);
}

template <typename T>
requires std::unsigned_integral<T>
T diffPosition(T a, T b) {
    auto differingBit = a ^ b;
    assert(std::popcount(differingBit) == 1);
    return std::countr_zero(differingBit);
}

#endif //SMS_GRAY_CODE_HPP
