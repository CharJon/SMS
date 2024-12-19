#ifndef GAWS24_CHRONO_IO_HPP
#define GAWS24_CHRONO_IO_HPP

#include <chrono>
#include <concepts>
#include <iostream>
#include <type_traits>

template <typename T>
    requires std::same_as<T, std::nano>
inline std::string suffix() {
    return "ns";
}

template <typename T>
    requires std::same_as<T, std::micro>
inline std::string suffix() {
    return "µs";
}

template <typename Ratio>
    requires std::same_as<Ratio, std::milli>
inline std::string suffix() {
    return "ms";
}

// no ratio means seconds
template <typename Rep>
std::stringstream &operator<<(std::stringstream &os, const std::chrono::duration<Rep> &duration) {
    os << duration.count() << "s";
    return os;
}

template <typename Rep>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep> &duration) {
    os << duration.count() << "s";
    return os;
}

// ratio allows to deduce unit
template <typename Rep, typename Ratio>
std::stringstream &operator<<(std::stringstream &os, const std::chrono::duration<Rep, Ratio> &duration) {
    os << duration.count() << suffix<Ratio>();
    return os;
}

template <typename Rep, typename Ratio>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep, Ratio> &duration) {
    os << duration.count() << suffix<Ratio>();
    return os;
}

namespace sms {

/**
 * A timer for running time measurements.
 * The timer is started on construction.
 */
class Timer {
public:
    using clock = std::chrono::steady_clock;

    Timer() { restart(); };

    clock::time_point restart();

    clock::time_point stop();

    clock::duration elapsed() const;

    clock::time_point startTime() const;

    clock::time_point stopTime() const;

    bool isRunning() const { return running_; }

protected:
    bool running_;
    clock::time_point started_;
    clock::time_point stopped_;

    clock::time_point stopTimeOrNow() const;
};
} // namespace sms

#endif // GAWS24_CHRONO_IO_HPP
