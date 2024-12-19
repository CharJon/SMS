#include "sms/auxiliary/cl_parser.hpp"

#include "cxxopts.hpp"
#include <chrono>
#include <optional>

#define INSTANCE "instance"
#define FILEFORMAT "filetype"
#define INSTTYPE "insttype"
#define SOLVER "solver"
#define SEED "seed"
#define TIMELIMIT "timelimit"
#define NODELIMIT "nodelimit"
#define PARAMFILE "paramfile"
#define STATSFILE "statsfile"
#define STATSSOLVERFILE "stats-solver"
#define LOGSOLVERFILE "log-solver"
#define HISTORYFILE "history"
#define SOLUTIONFILE "solution"
#define WARMSTARTSOLUTIONFILE "warm-start-solution"
#define THREADS "threads"
#define PRESOLVE "presolve-emphasis"
#define OLDPRESOLVE "presolve-old"
#define TRIANGLEEMPHASIS "triangle-emphasis"

cxxopts::Options CLParser::createParser() {
    cxxopts::Options options("SMS", "Max cut solver based on SCIP  ");

    options.add_options()
        // positional
        (INSTANCE, "Instance file name", cxxopts::value<std::string>()) //
        (FILEFORMAT, "File format (optional). Deduced from file ending if not present", cxxopts::value<std::string>())
        // pseudo-optional (have default values)
        (SOLVER, "Solver to use: (sms | none)", cxxopts::value<std::string>()->default_value("sms"))        //
        (SEED, "Random seed", cxxopts::value<int>()->default_value("0"))                                    //
        (TIMELIMIT, "Maximum time in seconds", cxxopts::value<int>()->default_value("2000000000"))          //
        (NODELIMIT, "Maximum number of nodes", cxxopts::value<int>()->default_value("-1"))                  //
        (PARAMFILE, "Parameter setting file to use", cxxopts::value<std::string>()->default_value(""))      //
        (PRESOLVE, "Presolve emphasis", cxxopts::value<int>()->default_value("1"))                          //
        (OLDPRESOLVE, "Only use old sota presolve", cxxopts::value<bool>())                                 //
        (TRIANGLEEMPHASIS, "Emphasis for triangles in presolve", cxxopts::value<int>()->default_value("3")) //
        // real optional
        (INSTTYPE, "Type of instance to interpret input. Possible values: bq, mc, fi, sg, gsg",
         cxxopts::value<std::string>())                                                            //
        (STATSFILE, "File to store compact stats in json format", cxxopts::value<std::string>())   //
        (STATSSOLVERFILE, "File to store solver specific stats in", cxxopts::value<std::string>()) //
        (LOGSOLVERFILE, "File to store solver log in", cxxopts::value<std::string>())              //
        (HISTORYFILE, "", cxxopts::value<std::string>())(SOLUTIONFILE, "Solution file",
                                                         cxxopts::value<std::string>())                //
        (WARMSTARTSOLUTIONFILE, "Solution file to use as a warm-start", cxxopts::value<std::string>()) //
        (THREADS, "Number of threads", cxxopts::value<int>()->default_value("1"))                      //
        ;
    options.parse_positional({INSTANCE, FILEFORMAT});

    return options;
}

CLParser::CLParser(int argc, char **argv) {
    parseResult_ = createParser().parse(argc, argv);
}

std::optional<std::string> CLParser::getFileName() const {
    if (parseResult_.count(INSTANCE) == 0) {
        return {};
    } else {
        return parseResult_[INSTANCE].as<std::string>();
    }
}

std::string CLParser::getFileFormat() const {
    if (parseResult_.count(FILEFORMAT) == 0) {
        uint64_t positionLastDot = getFileName().value().find_last_of(".");
        return getFileName().value().substr(positionLastDot + 1);
    } else {
        return parseResult_[FILEFORMAT].as<std::string>();
    }
}

std::string CLParser::getSolver() const {
    return parseResult_[SOLVER].as<std::string>();
}

int CLParser::getSeed() const {
    return parseResult_[SEED].as<int>();
}

std::chrono::seconds CLParser::getTimeLimit() const {
    return std::chrono::seconds(parseResult_[TIMELIMIT].as<int>());
}

std::optional<std::string> CLParser::getWarmStartSolutionFile() const {
    if (parseResult_.count(WARMSTARTSOLUTIONFILE) == 0)
        return std::nullopt;

    return parseResult_[WARMSTARTSOLUTIONFILE].as<std::string>();
}

std::optional<std::string> CLParser::getSolutionFile() const {
    if (parseResult_.count(SOLUTIONFILE) == 0)
        return std::nullopt;

    return parseResult_[SOLUTIONFILE].as<std::string>();
}

std::optional<std::string> CLParser::getHistoryFile() const {
    if (parseResult_.count(HISTORYFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[HISTORYFILE].as<std::string>();
}

std::optional<std::string> CLParser::getLogSolverFile() const {
    if (parseResult_.count(LOGSOLVERFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[LOGSOLVERFILE].as<std::string>();
}

std::optional<std::string> CLParser::getStatsSolverFile() const {
    if (parseResult_.count(STATSSOLVERFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[STATSSOLVERFILE].as<std::string>();
}

std::optional<std::string> CLParser::getStatsFile() const {
    if (parseResult_.count(STATSFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[STATSFILE].as<std::string>();
}

std::string CLParser::getParamFile() const {
    return parseResult_[PARAMFILE].as<std::string>();
}

int CLParser::getNodeLimit() const {
    return parseResult_[NODELIMIT].as<int>();
}

int CLParser::getThreads() const {
    return parseResult_[THREADS].as<int>();
}

std::string CLParser::getInstType() const {
    if (parseResult_.count(INSTTYPE) == 0) {
        // default should be mc, except if file type is bq or qplib
        if ((getFileFormat() == "bq") or (getFileFormat() == "qplib"))
            return "bq";
        else
            return "mc";
    } else {
        return parseResult_[INSTTYPE].as<std::string>();
    }
}

int CLParser::getPresolveEmphasis() const {
    return parseResult_[PRESOLVE].as<int>();
}

bool CLParser::getOldPresolve() const {
    return parseResult_[OLDPRESOLVE].as<bool>();
}

int CLParser::getTriangleEmphasis() const {
    return parseResult_[TRIANGLEEMPHASIS].as<int>();
}
