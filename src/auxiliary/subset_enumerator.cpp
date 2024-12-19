#include "sms/auxiliary/subset_enumerator.hpp"

#include <stdexcept>

namespace sms {

SubsetEnumerator::SubsetEnumerator(int64_t supersetSize, int64_t subsetSize)
    : supersetSize_(supersetSize), subsetSize_(subsetSize) {
    if (supersetSize_ > 63) { // first bit is needed for simple termination
        throw std::invalid_argument("Super set size too large for subset enumerator");
    }
    if (subsetSize_ > supersetSize_) {
        throw std::invalid_argument("Subset size larger than super set size");
    }

    currentSubset_ = firstSubset(); // all ones to the right
}

} // namespace sms