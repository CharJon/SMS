#ifndef SMS_SUBSET_ENUMERATOR_HPP
#define SMS_SUBSET_ENUMERATOR_HPP

#include <cassert>
#include <cstdint>

namespace sms {
/*
 * Encodes a subset as an 64 bit unsigned integer.
 * See here: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
 */
class SubsetEnumerator {

public:
    explicit SubsetEnumerator(int64_t supersetSize, int64_t subsetSize);

    // return true if all subsets got enumerated (and returned via calls to next)
    bool done() const {
        // done when all ones are to the left
        return currentSubset_ > (firstSubset() << (supersetSize_ - subsetSize_));
    }

    // return the next subset
    uint64_t next() {
        uint64_t old = currentSubset_;
        uint64_t t = (currentSubset_ | (currentSubset_ - 1)) + 1;
        currentSubset_ = t | ((((t & -t) / (currentSubset_ & -currentSubset_)) >> 1) - 1);
        return old;
    }

    // true if the current subset contains the given element
    bool contains(int elementPos) const {
        assert(!done());
        return currentSubset_ & (1 << elementPos);
    }

    uint64_t firstSubset() const { return (static_cast<uint64_t>(1) << subsetSize_) - 1; }

private:
    int64_t supersetSize_;
    int64_t subsetSize_;
    uint64_t currentSubset_;
};

} // namespace sms

#endif // SMS_SUBSET_ENUMERATOR_HPP
