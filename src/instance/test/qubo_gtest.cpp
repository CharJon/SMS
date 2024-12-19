#include <gtest/gtest.h>

#include "sms/instance/qubo.hpp"

TEST(QuBO, createClass1) {
    sms::QUBO qubo(5);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 5);
}

TEST(QuBO, createClass2) {
    sms::QUBO qubo(5);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 5);

    qubo.setValue(1, 1, 10.0);
    ASSERT_EQ(qubo.getValue(1, 1), 10.0);

    auto info = qubo.getInstanceInformation();

    ASSERT_EQ(info["dimension"].get<int>(), 5);
    ASSERT_EQ(info["scaling factor"].get<double>(), 1);
    ASSERT_EQ(info["offset"].get<double>(), 0);
}

TEST(QuBO, resetToZero1) {
    sms::QUBO qubo(100);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 100);

    qubo.setValue(1, 1, 7.3);
    ASSERT_EQ(qubo.getValue(1, 1), 7.3);

    qubo.resetToZero();
    for (int i = 0; i < qubo.getDim(); i++) {
        for (int j = 0; j < qubo.getDim(); j++) {
            ASSERT_EQ(qubo.getValue(i, j), 0.0);
        }
    }
}

TEST(QuBO, quboEval1) {
    sms::QUBO qubo(100);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 100);

    qubo.resetToZero();

    double solVector[100];
    std::fill_n(solVector, 100, 1.0);

    ASSERT_EQ(qubo.getSolutionValue(solVector), 0.0);
}

TEST(QuBO, quboEvalBLAS1) {
    sms::QUBO qubo(100);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 100);

    qubo.resetToZero();

    double solVector[100];
    std::fill_n(solVector, 100, 1.0);

    ASSERT_EQ(qubo.getSolutionValueBLAS(solVector), 0.0);
}

TEST(QuBO, Copy) {
    sms::QUBO qubo(100);
    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 100);

    qubo.setValue(1, 1, 7.3);
    ASSERT_EQ(qubo.getValue(1, 1), 7.3);

    sms::QUBO qubo2(qubo);
    ASSERT_TRUE(qubo2.isValid());
    ASSERT_EQ(qubo2.getDim(), 100);

    ASSERT_EQ(qubo2.getValue(1, 1), 7.3);

    qubo2.setValue(1, 1, 8.3);

    ASSERT_EQ(qubo2.getValue(1, 1), 8.3);
    ASSERT_EQ(qubo.getValue(1, 1), 7.3);
}

TEST(QuBO, fromBqFile) {
    auto qubo = sms::bqFileToQubo("test/data/square.bq");

    ASSERT_TRUE(qubo.isValid());
    ASSERT_EQ(qubo.getDim(), 4);
    ASSERT_EQ(qubo.getValue(0, 1), 0.5);
    ASSERT_EQ(qubo.getValue(0, 3), -1.2);
    ASSERT_EQ(qubo.getValue(1, 2), -2);
    ASSERT_EQ(qubo.getValue(2, 3), 0.6);
    ASSERT_EQ(qubo.getValue(3, 3), 0.5);
}