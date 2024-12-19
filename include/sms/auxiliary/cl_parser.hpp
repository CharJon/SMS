#ifndef SMS_CL_PARSER_HPP
#define SMS_CL_PARSER_HPP

#include "cxxopts.hpp"
#include <chrono>
#include <optional>

class CLParser {
public:
    CLParser(int argc, char **argv);

    std::optional<std::string> getFileName() const;

    std::string getFileFormat() const;

    std::string getInstType() const;

    std::string getSolver() const;

    int getSeed() const;

    std::chrono::seconds getTimeLimit() const;

    int getNodeLimit() const;

    std::string getParamFile() const;

    std::optional<std::string> getStatsFile() const;

    std::optional<std::string> getStatsSolverFile() const;

    std::optional<std::string> getLogSolverFile() const;

    std::optional<std::string> getHistoryFile() const;

    std::optional<std::string> getSolutionFile() const;

    std::optional<std::string> getWarmStartSolutionFile() const;

    int getThreads() const;

    int getPresolveEmphasis() const;

    bool getOldPresolve() const;

    int getTriangleEmphasis() const;

private:
    cxxopts::ParseResult parseResult_;

    static cxxopts::Options createParser();
};

#endif // SMS_CL_PARSER_HPP
