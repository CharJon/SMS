#ifndef SMS_ACTIVE_ELEMENTS_HPP
#define SMS_ACTIVE_ELEMENTS_HPP

#include <cassert>
#include <concepts>
#include <limits>
#include <vector>

namespace sms {

/**
 * @brief Keeps track of n integer elements between [0, n-1].
 * Elements can be activated in O(1) and deactivated in O(1), and their status can be checked in O(1).
 * Also one active element can be popped in O(1).
 * @tparam T The type of the elements.
 */
template <std::unsigned_integral T>
class ActiveElementsStack {
public:
    explicit ActiveElementsStack(T n) : n_(n), stack_(), posInStack_(n_, n_) {
        assert(n_ < std::numeric_limits<T>::max()); // We need to reserve the max to flag inactive elements
        stack_.reserve(n_);
    }

    bool isActive(T x) {
        assert(x < n_);
        return posInStack_[x] != n_;
    }

    void activate(T x) {
        assert(x < n_);
        if (!isActive(x)) {
            posInStack_[x] = stack_.size();
            stack_.push_back(x);
        }
    }

    void deactivate(T x) {
        assert(x < n_);
        if (isActive(x)) {
            auto back = stack_.back();
            posInStack_[back] = posInStack_[x];
            stack_[posInStack_[back]] = back;
            stack_.pop_back();
            posInStack_[x] = n_;
        }
    }

    T popBack() {
        assert(!empty());
        auto last = stack_.back();
        posInStack_[last] = n_;
        stack_.pop_back();
        return last;
    }

    bool empty() { return stack_.empty(); }

private:
    T n_;
    std::vector<T> stack_;
    std::vector<T> posInStack_;
};

} // namespace sms

#endif // SMS_ACTIVE_ELEMENTS_HPP
