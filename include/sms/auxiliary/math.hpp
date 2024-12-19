#ifndef SMS_MATH_HPP
#define SMS_MATH_HPP

#include "cblas.h"
#include <array>
#include <climits>
#include <cmath>
#include <numeric>

/**
 * @file math.hpp
 * @brief Math utilities
 */

namespace sms {

// Divide an integer by two and ceil the result
template <std::integral T>
inline T divTwoCeil(T a) {
    return (a + T{1}) / 2;
}

// All positive integers equal to and smaller than this value can be exactly represented as T
// in  iec 559 (same as IEEE 754) floats 32 or 64 bit
template <std::floating_point T>
constexpr int64_t maxSafeInteger() {
    static_assert(sizeof(T) * CHAR_BIT == 32 || sizeof(T) * CHAR_BIT == 64);
    static_assert(std::numeric_limits<T>::is_iec559);
    constexpr int64_t val = 1;
    if (sizeof(T) * CHAR_BIT == 64)
        return val << 53;
    else if (sizeof(T) * CHAR_BIT == 32)
        return val << 23;
}

/**
 * @brief Computes alpha * x^T * Q * x
 * @param Q Q matrix in row major order
 * @param x x vector
 * @param scalingFactor alpha
 * @return x^T * Q * x
 */
inline double quboEvaluation(const double *Q, const double *x, int dim, double scalingFactor) {
    auto *y = new double[dim];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, scalingFactor, Q, dim, x, 1, 0.0, y, 1);
    double res = cblas_ddot(dim, x, 1, y, 1);
    delete[] y;
    return res;
}

/**
 * @brief Computes x^T * Q * x
 * @param Q Q matrix in row major order
 * @param x x vector
 * @return x^T * Q * x
 */
inline double quboEvaluation(const double *Q, const double *x, int dim) {
    return quboEvaluation(Q, x, dim, 1.0);
}

/**
 * @brief Checks if two doubles are equal within a given absolute tolerance
 * @param a first double
 * @param b second double
 * @param eps epsilon
 * @return true if a and b are equal within eps, false otherwise
 */
inline bool equalAbsEps(double a, double b, double eps) {
    return std::fabs(a - b) <= eps;
}

inline bool equalRelEps(double a, double b, double eps) {
    return std::fabs(a - b) <= eps * std::max(std::fabs(a), std::fabs(b));
}

inline bool equalEps(double a, double b, double absEps, double relEps) {
    auto lhs = std::fabs(a - b);
    auto rhs = std::max(absEps, relEps * std::max(std::fabs(a), std::fabs(b)));
    return lhs <= rhs;
}

inline bool equalEps(double a, double b, double eps) {
    return equalEps(a, b, eps, eps);
}

/**
 * @brief Checks if a double is less than (<) another within a given absolute tolerance
 * @param a first double
 * @param b second double
 * @param eps epsilon
 * @return true if a < b within eps, false otherwise
 */
inline bool ltEps(double a, double b, double eps) {
    return a - b < -eps;
}

/**
 * @brief Checks if a double is less than or equal to (<=) another within a given absolute tolerance
 * @param a first double
 * @param b second double
 * @param eps epsilon
 * @return true if a <= b within eps, false otherwise
 */
inline bool leEps(double a, double b, double eps) {
    return a - b <= eps;
}

/**
 * @brief Checks if a double is greater than (>) another within a given absolute tolerance
 * @param a first double
 * @param b second double
 * @param eps epsilon
 * @return true if a > b within eps, false otherwise
 */
inline bool gtEps(double a, double b, double eps) {
    return a - b > eps;
}

/**
 * @brief Checks if a double is greater than or equal to (>=) another within a given absolute tolerance
 * @param a first double
 * @param b second double
 * @param eps epsilon
 * @return true if a >= b within eps, false otherwise
 */
inline bool geEps(double a, double b, double eps) {
    return a - b >= -eps;
}

/**
 * @brief Checks if a double is zero within a given absolute tolerance
 * @param a double to check
 * @param eps epsilon
 * @return true if a is zero within eps, false otherwise
 */
inline bool zeroEps(double a, double eps) {
    return std::fabs(a) < eps;
}

// Fast way to check if a floating point number is an integer. Infinity returns true!
template <std::floating_point T>
bool isInteger(T a) {
    return std::rint(a) == a;
}

// Fast way to check if a floating point number is epsilon close to an integer. Infinity returns true!
template <std::floating_point T>
bool isNearlyInteger(T a, T eps = 1e-6) {
    T intPart;
    return zeroEps(std::modf(a, &intPart), eps);
}

/*
 * Discussion from stackoverflow:
 * https://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
 */
constexpr double fastCos(double phi) {
    constexpr double tp = 1. / (2. * std::numbers::pi);
    phi *= tp;
    phi -= .25 + std::floor(phi + .25);
    phi *= 16. * (std::abs(phi) - .5);
    phi += .225 * phi * (std::abs(phi) - 1.);
    return phi;
}

constexpr double fastSin(const double phi) {
    return fastCos(phi - std::numbers::pi / 2.);
}

// Tries to find a rational representation of @val with a denominator smaller than @maxDenominator with an error of @eps
// Returns a pair of the numerator and denominator of the rational representation
inline std::pair<int, int> rationalRepresentation(double val, int maxDenominator, double eps = 10e-6) {
    if (maxDenominator < 1)
        throw std::runtime_error("maxDominator should be at least 1");
    if (isNearlyInteger(val, eps))
        return {static_cast<int>(val), 1};

    int sign = val >= 0 ? 1 : -1;
    val = std::abs(val);

    constexpr std::array<int, 19> simpleDenoms{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 25};

    // try the simple denominators first
    int scaling = 1;

    while (scaling < maxDenominator) {
        for (auto d : simpleDenoms) {
            double denom = scaling * d;
            if (denom <= maxDenominator) {
                double num = std::round(val * denom);
                if (equalEps(val, num / denom, eps)) {
                    auto gcd = std::gcd(static_cast<int>(num), static_cast<int>(denom));
                    return {sign * (num / gcd), denom / gcd};
                }
            }
        }
        scaling *= 10;
    }

    // approximate the rational representation using continued fractions,
    // see https://en.wikipedia.org/wiki/Continued_fraction#Calculating_continued_fraction_representations
    double x = val;
    int floorX = std::floor(x); // integer part of val

    int num0 = floorX;
    int num1 = 1;
    int denom0 = 1;
    int denom1 = 0;

    while (denom0 < maxDenominator) {
        x = 1.0 / (x - floorX); // inverse of (remaining) fractional part
        floorX = std::floor(x);
        int curNum = num0;
        int curDenom = denom0;

        num0 = floorX * num0 + num1;
        denom0 = floorX * denom0 + denom1;
        // num0 / denom0 = continued fraction
        if (equalEps(val, static_cast<double>(num0) / static_cast<double>(denom0), eps))
            break;
        if (equalEps(val, (num0 - 1.0) / static_cast<double>(denom0), eps)) {
            num0 = num0 - 1;
            break;
        }
        if (equalEps(val, (num0 + 1.0) / static_cast<double>(denom0), eps)) {
            num0 = num0 + 1;
            break;
        }

        num1 = curNum;
        denom1 = curDenom;
    }

    if (denom0 <= maxDenominator) {
        auto gcd = std::gcd(num0, denom0);
        return {sign * (num0 / gcd), denom0 / gcd};
    }
    return {0, -1};
}

} // namespace sms

#endif // SMS_MATH_HPP
