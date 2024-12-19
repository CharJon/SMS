#ifndef SMS_SG_PARSER_HPP
#define SMS_SG_PARSER_HPP

#include <optional>
#include <string>

#include "sms/instance/gridded_ising.hpp"
#include "sms/io/io.hpp"

namespace sms {

class SgParser {
public:
    explicit SgParser(const std::string &path, double scaling = 1.0, double offset = 0.0);

    Ising getInstance();

private:
    NetworKit::Graph graph_;
    double scaling_;
    double offset_;
};

class GsgParser {
public:
    explicit GsgParser(const std::string &path, double scaling = 1.0, double offset = 0.0);

    std::optional<torusIsing> getTorus() const;

    std::optional<gridIsing> getGrid() const;

    bool isTorus() const;

private:
    NetworKit::Graph graph_;
    double scaling_;
    double offset_;
    int dim_;
    int gridLength_;
    bool torus_;

    void isValid();
};

} // namespace sms

#endif // SMS_SG_PARSER_HPP
