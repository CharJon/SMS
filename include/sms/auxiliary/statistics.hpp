#ifndef SMS_STATISTICS_HPP
#define SMS_STATISTICS_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include "nlohmann/json.hpp"

#include "networkit/graph/Graph.hpp"

namespace sms {

template <typename T>
    requires std::integral<T> || std::floating_point<T>
T getMedian(std::vector<T> &values) {
    size_t n = values.size();
    size_t medianPos = (n - 1) / 2;
    std::nth_element(values.begin(), values.begin() + medianPos, values.end());
    return values[medianPos];
}

template <typename T>
    requires std::integral<T> || std::floating_point<T>
double getMean(std::vector<T> &values) {
    auto n = static_cast<double>(values.size());
    return std::accumulate(values.begin(), values.end(), 0.0, [&n](double cma, double i) { return cma + (i / n); });
}

template <typename T>
    requires std::integral<T> || std::floating_point<T>
double getMean2(std::vector<T> &values) {
    return std::accumulate(values.begin(), values.end(), 0.0,
                           [n = 0](auto cma, auto i) mutable { return cma + (i - cma) / ++n; });
}

template <typename T>
    requires std::integral<T> || std::floating_point<T>
double getStdDev(std::vector<T> &values, double mean) {
    double stdDev = 0;
    for (T x : values) {
        stdDev += (x - mean) * (x - mean);
    }
    return std::sqrt(stdDev / values.size());
}

template <std::integral T>
T getMode(std::vector<T> &values) {
    assert(not values.empty());

    std::map<T, int> modeMap;
    for (int x : values) {
        modeMap[x] = 0;
    }

    int mode = 0;
    int maxModeCount = 0;

    for (int x : values) {
        modeMap[x] += 1;
        if (modeMap[x] > maxModeCount) {
            mode = x;
            maxModeCount = modeMap[x];
        }
    }

    return mode;
}

template <std::integral T>
class IntStatistics {
public:
    explicit IntStatistics(std::vector<T> values) {
        assert(not values.empty());

        min = *std::min_element(values.begin(), values.end());
        max = *std::max_element(values.begin(), values.end());
        mean = getMean(values);
        median = getMedian(values);
        stdDev = getStdDev(values, mean);
        mode = getMode(values);
    }

    nlohmann::ordered_json asJson() {
        return nlohmann::ordered_json{
            {"min", min}, {"max", max}, {"mean", mean}, {"median", median}, {"standard_dev", stdDev}, {"mode", mode}};
    }

    T min;
    T max;
    double mean;
    T median;
    double stdDev;
    T mode;
};

template <std::floating_point T>
class FloatStatistics {
public:
    explicit FloatStatistics(std::vector<T> values) {
        assert(not values.empty());

        min = *std::min_element(values.begin(), values.end());
        max = *std::max_element(values.begin(), values.end());
        mean = getMean(values);
        median = getMedian(values);
        stdDev = getStdDev(values, mean);
    }

    nlohmann::ordered_json asJson() {
        return nlohmann::ordered_json{
            {"min", min}, {"max", max}, {"mean", mean}, {"median", median}, {"standard_dev", stdDev}};
    }

    T min;
    T max;
    T mean;
    T median;
    T stdDev;
};

} // namespace sms

#endif // SMS_STATISTICS_HPP