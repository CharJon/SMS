#include "cxxopts.hpp"
#include <iostream>

#include "sms/io/dimacs_reader.hpp"
#include "sms/io/io.hpp"
#include "sms/io/ising_reader.hpp"
#include "sms/io/qplib_parser.hpp"

cxxopts::Options cliParser() {
    cxxopts::Options options("SMS", "Max cut solver based on SCIP  ");

    options.add_options()
        // positional
        ("file", "Instance file name", cxxopts::value<std::string>())(
            "type", "Type (optional)", cxxopts::value<std::string>())("hidepaths", "Hide full paths in error logs",
                                                                      cxxopts::value<bool>()->default_value("false"));

    options.parse_positional({"file", "type"});

    return options;
}

int main(int args, char *argv[]) {

    auto opts = cliParser();

    auto result = opts.parse(args, argv);

    if (result.count("file") == 0) {
        std::string usage = "Usage:\n"
                            "1 - path to file to check\n"
                            "2 - (optional) type of check to do, either mc or bq\n"
                            "3 - (optional) hide full path of file in error log\n";
        throw std::runtime_error(usage);
    }

    std::string type;
    std::string file = result["file"].as<std::string>();

    if (result.count("type") == 0) {
        uint64_t position_last_dot = file.find_last_of('.');
        type = file.substr(position_last_dot + 1);
    } else {
        type = result["type"].as<std::string>();
    }

    try {
        if (type == "mc" or type == "sg") {
            fileToMcObj(file);
            return EXIT_SUCCESS;

        } else if (type == "bq") {
            fileToBqObj(file);
            return EXIT_SUCCESS;

        } else if (type == "gsg") {
            auto parser = sms::GsgParser(file);
            return EXIT_SUCCESS;

        } else if (type == "dimacs") {
            auto parser = sms::DimacsReader(file);
            return EXIT_SUCCESS;

        } else if (type == "qplib") {
            auto parser = sms::QPLibParser(file);
            return EXIT_SUCCESS;

        } else {
            throw std::runtime_error("Unknown check type: " + type);
        }
    } catch (const std::exception &e) {
        std::string filename = file;

        if (result["hidepaths"].count() == 1) {
            uint64_t position_last_slash = file.find_last_of('/');
            filename = file.substr(position_last_slash + 1);
        }

        std::cerr << "Error with file " << filename << std::endl;
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}