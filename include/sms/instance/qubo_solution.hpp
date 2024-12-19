#ifndef SMS_QUBO_SOLUTION_HPP
#define SMS_QUBO_SOLUTION_HPP

#include <cassert>
#include <fstream>
#include <vector>

#include "nlohmann/json.hpp"

#include "sms/auxiliary/bipartition.hpp"
#include "sms/instance/qubo.hpp"

namespace sms {

using var_t = uint64_t;
using assignment_t = uint8_t;

class QuboSolution {
public:
    explicit QuboSolution(const QUBO *qubo);

    // not implemented yet
    explicit QuboSolution(const QUBO *qubo, const std::string &file);

    // not implemented yet
    explicit QuboSolution(const QUBO *qubo, const std::vector<assignment_t> &assignments);

    void assignTo0(var_t u);

    void assignTo1(var_t u);

    void allVarsTo0();

    void allVarsTo1();

    double getQuboValue() const;

    std::vector<var_t> getAssignment0() const;

    std::vector<var_t> getAssignment1() const;

    assignment_t getAssignment(var_t u) const;

    void flipAssignment();

    void saveToFile(const std::string &) const;

    bool isValid() const;

    const uint8_t kUNASSIGNED = 2; // used for unassigned nodes

private:
    const QUBO *const qubo_;
    std::vector<assignment_t> assignment_;
};

} // namespace sms

#endif // SMS_QUBO_SOLUTION_HPP
