#include <random>

#include <gtest/gtest.h>

#include "sms/auxiliary/math.hpp"

TEST(Math, isInteger) {
    EXPECT_TRUE(sms::isInteger(1.0));
    EXPECT_TRUE(sms::isInteger(2.0));
    EXPECT_TRUE(sms::isInteger(2.0 * 5.0));
    EXPECT_FALSE(sms::isInteger(1.00001));
    EXPECT_FALSE(sms::isInteger(1.0 / 3.0));
    EXPECT_FALSE(sms::isInteger(std::numeric_limits<double>::quiet_NaN()));
    EXPECT_TRUE(sms::isInteger(std::numeric_limits<double>::infinity()));
    EXPECT_FALSE(sms::isInteger(std::numeric_limits<double>::signaling_NaN()));
    EXPECT_FALSE(sms::isInteger(std::numeric_limits<double>::denorm_min()));
}

TEST(Math, isNearlyInteger) {
    EXPECT_TRUE(sms::isNearlyInteger(1.0));
    EXPECT_TRUE(sms::isNearlyInteger(2.0));
    EXPECT_TRUE(sms::isNearlyInteger(2.0 * 5.0));

    EXPECT_FALSE(sms::isNearlyInteger(1.00001));
    EXPECT_FALSE(sms::isNearlyInteger(1.0 / 3.0));

    double eps = 1e-6;
    EXPECT_TRUE(sms::isNearlyInteger(2.0 + eps, 2 * eps));

    EXPECT_FALSE(sms::isNearlyInteger(std::numeric_limits<double>::quiet_NaN()));
    EXPECT_TRUE(sms::isNearlyInteger(std::numeric_limits<double>::infinity()));
    EXPECT_FALSE(sms::isNearlyInteger(std::numeric_limits<double>::signaling_NaN()));
    EXPECT_TRUE(sms::isNearlyInteger(std::numeric_limits<double>::denorm_min()));
}

TEST(Math, BlasVectorMultMatrix) {
    double matA[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double x[3] = {1, 1, 1};
    double y[3] = {0, 0, 0};

    // Perform matrix-vector multiplication
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, &matA[0][0], 3, x, 1, 0.0, y, 1);

    ASSERT_EQ(y[0], 6);
    ASSERT_EQ(y[1], 15);
    ASSERT_EQ(y[2], 24);
}

TEST(Math, BlasQUBO) {
    double matA[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double x[3] = {1, 0, 1};
    double y[3] = {0, 0, 0};

    // Perform matrix-vector multiplication
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, &matA[0][0], 3, x, 1, 0.0, y, 1);
    auto res = cblas_ddot(3, x, 1, y, 1);

    ASSERT_EQ(res, 20);
}

TEST(Math, BlasQUBO_2) {
    double matA[2][2] = {{-1, -1}, {-1, -1}};
    double x[3] = {1, 1};
    double y[3] = {0, 0};

    // Perform matrix-vector multiplication
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2, 2, 1.0, &matA[0][0], 2, x, 1, 0.0, y, 1);
    auto res = cblas_ddot(2, x, 1, y, 1);

    ASSERT_EQ(res, -4);
}

TEST(Math, Numerics) {
    EXPECT_TRUE(sms::equalAbsEps(1.0, 1.0, 0.0001));
    EXPECT_TRUE(sms::equalAbsEps(1.0, 1.1, 0.10001));
    EXPECT_FALSE(sms::equalAbsEps(1.0, 1.1, 0.01));

    EXPECT_TRUE(sms::ltEps(1.0, 1.1, 0.001));
    EXPECT_FALSE(sms::ltEps(1.0, 1.09, 0.099));
    EXPECT_TRUE(sms::leEps(1.0, 1.1, 0.1));

    EXPECT_FALSE(sms::gtEps(1.0, 1.1, 0.1));
    EXPECT_TRUE(sms::gtEps(1.1, 1.04, 0.05));
    EXPECT_TRUE(sms::geEps(1.1, 1.19, 0.1));

    EXPECT_TRUE(sms::zeroEps(0.01, 0.1));

    EXPECT_FALSE(sms::equalAbsEps(100.0, 101.0, 0.1));
    EXPECT_TRUE(sms::equalRelEps(100.0, 101.0, 0.1));
    EXPECT_TRUE(sms::equalEps(100.0, 101.0, 0.01));
    EXPECT_FALSE(sms::equalEps(100.0, 101.0, 0.009));
}

TEST(Math, Limits) {
    EXPECT_EQ(sms::maxSafeInteger<float>(), pow(2, 23));
    EXPECT_EQ(sms::maxSafeInteger<double>(), pow(2, 53));

    auto doubleLimit = sms::maxSafeInteger<double>();
    EXPECT_EQ(doubleLimit, static_cast<int64_t>(static_cast<double>(doubleLimit)));
    EXPECT_NE(doubleLimit + 1, static_cast<int64_t>(static_cast<double>(doubleLimit + 1)));
}

TEST(Math, FastCos) {
    for (int i = 0; i < 100; i++) {
        std::cout << i * 0.1 << std::endl;
        double x = i * 0.1;

        double cosX = std::cos(x);
        double fastCosX = sms::fastCos(x);
        EXPECT_NEAR(cosX, fastCosX, 1e-2);

        double sinX = std::sin(x);
        double fastSinX = sms::fastSin(x);
        EXPECT_NEAR(sinX, fastSinX, 1e-2);
    }

    double pi = std::numbers::pi;
    std::array<double, 7> testAngles = {0, pi / 2, pi, 3 * pi / 2, 2 * pi, -pi / 2, -pi};
    for (double angle : testAngles) {
        std::cout << angle << std::endl;

        double cosX = std::cos(angle);
        double fastCosX = sms::fastCos(angle);
        EXPECT_NEAR(cosX, fastCosX, 1e-2);

        double sinX = std::sin(angle);
        double fastSinX = sms::fastSin(angle);
        EXPECT_NEAR(sinX, fastSinX, 1e-2);
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1000 * pi, 1000 * pi);

    for (int i = 0; i < 1000; i++) {
        double x = distribution(generator);
        std::cout << x << std::endl;

        double cosX = std::cos(x);
        double fastCosX = sms::fastCos(x);
        EXPECT_NEAR(cosX, fastCosX, 1e-2);

        double sinX = std::sin(x);
        double fastSinX = sms::fastSin(x);
        EXPECT_NEAR(sinX, fastSinX, 1e-2);
    }
}