#include "sms/auxiliary/chrono_io.hpp"

namespace sms {

Timer::clock::time_point Timer::restart() {
    started_ = clock::now();
    running_ = true;
    return started_;
}

Timer::clock::time_point Timer::stop() {
    stopped_ = clock::now();
    running_ = false;
    return stopped_;
}

Timer::clock::duration Timer::elapsed() const {
    return stopTimeOrNow() - started_;
}

Timer::clock::time_point Timer::startTime() const {
    return started_;
}

Timer::clock::time_point Timer::stopTime() const {
    return stopped_;
}

Timer::clock::time_point Timer::stopTimeOrNow() const {
    return running_ ? clock::now() : stopped_;
}

} // namespace sms