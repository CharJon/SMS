#include <filesystem>

#include <gtest/gtest.h>

#include "sms/io/ising_reader.hpp"

TEST(SGParser, SimpleParsing) {
    sms::SgParser parser("test/data/square.sg");

    auto inst = parser.getInstance();
    ASSERT_EQ(inst.getNumberOfSpins(), 4);
    ASSERT_EQ(inst.getNumberOfInteractions(), 4);
    ASSERT_EQ(inst.getSolutionValue({1, 1, 1, 1}), 7);
    ASSERT_EQ(inst.getSolutionValue({1, -1, 1, 1}), -3);
    ASSERT_EQ(inst.getSolutionValue({1, -1, -1, 1}), 25);
}

TEST(SGParser, PropblemticInstance) {
    ASSERT_THROW(sms::SgParser parser("test/data/mcsparse_online_instances/sg_error/isg_pm_81_9.sg"),
                 std::runtime_error);
}

TEST(GSGParser, SimpleParsing) {
    sms::GsgParser parser("test/data/simpleGrid.gsg");

    auto inst = parser.getTorus().value();
    ASSERT_EQ(inst.getNumberOfSpins(), 9);
    ASSERT_EQ(inst.getNumberOfInteractions(), 18);
    ASSERT_EQ(inst.getDim(), 2);
    ASSERT_EQ(inst.getGridLength(), 3);
    ASSERT_TRUE(inst.isConsistent());
}

TEST(GSGParser, SimpleParsingNonPeriodic) {
    sms::GsgParser parser("test/data/simpleGrid_non_periodic.gsg");

    auto inst = parser.getGrid().value();
    ASSERT_EQ(inst.getNumberOfSpins(), 9);
    ASSERT_EQ(inst.getNumberOfInteractions(), 12);
    ASSERT_EQ(inst.getDim(), 2);
    ASSERT_EQ(inst.getGridLength(), 3);
    ASSERT_TRUE(inst.isConsistent());
}

TEST(GSGParser, 3DParsing) {
    sms::GsgParser parser("test/data/t3g7_6666.gsg");

    auto inst = parser.getTorus().value();
    ASSERT_EQ(inst.getNumberOfSpins(), 343);
    ASSERT_EQ(inst.getNumberOfInteractions(), 1029);
    ASSERT_EQ(inst.getDim(), 3);
    ASSERT_EQ(inst.getGridLength(), 7);
    ASSERT_TRUE(inst.isConsistent());
}

TEST(GSGParser, WrongGrid) {
    ASSERT_ANY_THROW(sms::GsgParser parser("test/data/error_grid_undetermined.gsg"));
}

TEST(GSGParser, MissingEdge) {
    ASSERT_ANY_THROW(sms::GsgParser parser("test/data/error_missing_edge.gsg"));
}

TEST(GSGParser, TooFewEdges) {
    ASSERT_ANY_THROW(sms::GsgParser parser("test/data/error_wrong_number_of_edges.gsg"));
}

TEST(GSGParser, TooManyEdges) {
    ASSERT_ANY_THROW(sms::GsgParser parser("test/data/error_wrong_number_of_edges_2.gsg"));
}

TEST(gsgRead, McSparseOnlineError) {
    std::string path = "test/data/mcsparse_online_instances/gsg_error";
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".gsg") {
            ASSERT_ANY_THROW(sms::GsgParser parser(entry.path()));
        }
    }
}

TEST(gsgRead, McSparseOnlineValid) {
    std::string path = "test/data/mcsparse_online_instances/gsg_valid";
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".gsg") {
            ASSERT_NO_THROW(sms::GsgParser parser(entry.path()));
        }
    }
}

TEST(sgRead, McSparseOnlineError) {
    std::string path = "test/data/mcsparse_online_instances/sg_error";
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".sg") {
            ASSERT_ANY_THROW(sms::SgParser parser(entry.path()));
        }
    }
}

TEST(sgRead, McSparseOnlineValid) {
    std::string path = "test/data/mcsparse_online_instances/sg_valid";
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".sg") {
            ASSERT_NO_THROW(sms::SgParser parser(entry.path()));
        }
    }
}