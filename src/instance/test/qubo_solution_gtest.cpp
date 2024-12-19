#include <gtest/gtest.h>

#include "sms/instance/qubo.hpp"
#include "sms/instance/qubo_solution.hpp"

TEST(QuboSolutionClass, Create) {
    auto inst = sms::QUBO(4);

    sms::QuboSolution solution(&inst);
}

TEST(QuboSolutionClass, Init) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -3);

    sms::QuboSolution solution(&inst);

    solution.allVarsTo0();
    ASSERT_EQ(solution.getQuboValue(), 0);

    solution.allVarsTo1();
    ASSERT_EQ(solution.getQuboValue(), 1);
}

TEST(QuboSolutionClass, SquareUnweightedValue1) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst);

    solution.allVarsTo0();
    solution.assignTo1(2);
    ASSERT_EQ(solution.getQuboValue(), 1);
}

TEST(QuboSolutionClass, SquareWeightedCut1) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst);

    solution.allVarsTo0();
    solution.assignTo1(0);
    solution.assignTo1(2);
    ASSERT_EQ(solution.getQuboValue(), 2);
}

TEST(QuboSolutionClass, SquareFlip1) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst);

    solution.allVarsTo0();

    auto part0nodesBefore = solution.getAssignment0();
    auto part1nodesBefore = solution.getAssignment1();
    ASSERT_EQ(part0nodesBefore.size(), 4);
    ASSERT_EQ(part1nodesBefore.size(), 0);

    solution.flipAssignment();
    auto part0nodesAfter = solution.getAssignment0();
    auto part1nodesAfter = solution.getAssignment1();
    ASSERT_EQ(part0nodesAfter.size(), 0);
    ASSERT_EQ(part1nodesAfter.size(), 4);
}

TEST(QuboSolutionClass, LineTest) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst);

    solution.allVarsTo0();
    solution.assignTo1(2);
    solution.assignTo1(3);
    ASSERT_EQ(solution.getQuboValue(), 2);
}

TEST(QuboSolutionClass, LoadTest) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst, "test/data/qubo_solution_test.json");

    ASSERT_EQ(solution.getQuboValue(), 2);
    ASSERT_TRUE(solution.isValid());
}

TEST(QuboSolutionClass, LoadTestEmptyAssignnment0) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst, "test/data/qubo_solution_test_empty_0.json");

    ASSERT_EQ(solution.getQuboValue(), 0);
    ASSERT_TRUE(solution.isValid());
}

TEST(QuboSolutionClass, LoadTestEmptyAssignnment1) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst, "test/data/qubo_solution_test_empty_1.json");

    ASSERT_EQ(solution.getQuboValue(), 0);
    ASSERT_TRUE(solution.isValid());
}

TEST(QuboSolutionClass, LoadTestEmptyBoth) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    sms::QuboSolution solution(&inst, "test/data/qubo_solution_test_double_empty.json");

    ASSERT_FALSE(solution.isValid());
}

TEST(QuboSolutionClass, LoadTestDoubleDef) {
    auto inst = sms::QUBO(4);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);
    inst.setValue(3, 3, 1);
    inst.setValue(3, 1, -4);

    ASSERT_DEBUG_DEATH(sms::QuboSolution solution(&inst, "test/data/qubo_solution_test_double_definition.json"), "");
}

TEST(QuboSolutionClass, LoadTestOutOfBounds) {
    auto inst = sms::QUBO(3);
    inst.setValue(0, 0, 1);
    inst.setValue(1, 1, 1);
    inst.setValue(2, 2, 1);

    ASSERT_DEBUG_DEATH(sms::QuboSolution solution(&inst, "test/data/qubo_solution_test.json"), "");
}