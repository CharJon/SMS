#include "sms/instance/qubo_solution.hpp"

#include <cassert>
#include <fstream>
#include <vector>

#include "nlohmann/json.hpp"

namespace sms {

#define ASSIGN_0_NAME "assignment_0"
#define ASSIGN_1_NAME "assignment_1"

QuboSolution::QuboSolution(const QUBO *qubo) : qubo_(qubo) {
    assignment_ = std::vector<assignment_t>(qubo_->getDim(), kUNASSIGNED);
}

void QuboSolution::assignTo0(var_t u) {
    assert(qubo_->getDim() > u);
    assignment_[u] = 0;
}

void QuboSolution::assignTo1(var_t u) {
    assert(qubo_->getDim() > u);
    assignment_[u] = 1;
}

void QuboSolution::allVarsTo0() {
    for (int i = 0; i < qubo_->getDim(); ++i) {
        assignment_[i] = 0;
    }
}

void QuboSolution::allVarsTo1() {
    for (int i = 0; i < qubo_->getDim(); ++i) {
        assignment_[i] = 1;
    }
}

double QuboSolution::getQuboValue() const {
    auto solArray = new double[qubo_->getDim()];

    for (int i = 0; i < qubo_->getDim(); ++i) {
        solArray[i] = assignment_[i];
    }

    return qubo_->getSolutionValue(solArray);
}

std::vector<var_t> QuboSolution::getAssignment0() const {
    std::vector<var_t> assign0vars{};
    for (int u = 0; u < qubo_->getDim(); ++u) {
        if (assignment_[u] == 0)
            assign0vars.push_back(u);
    }
    return assign0vars;
}

std::vector<var_t> QuboSolution::getAssignment1() const {
    std::vector<var_t> assign1vars{};
    for (int u = 0; u < qubo_->getDim(); ++u) {
        if (assignment_[u] == 1)
            assign1vars.push_back(u);
    }
    return assign1vars;
}

assignment_t QuboSolution::getAssignment(var_t u) const {
    assert(qubo_->getDim() > u);
    return assignment_[u];
}

void QuboSolution::flipAssignment() {
    for (int i = 0; i < qubo_->getDim(); ++i) {
        assignment_[i] = assignment_[i] == kUNASSIGNED ? kUNASSIGNED : 1 - assignment_[i];
    }
}

void QuboSolution::saveToFile(const std::string &path) const {
    nlohmann::json save;

    save[ASSIGN_0_NAME] = getAssignment0();
    save[ASSIGN_1_NAME] = getAssignment1();

    std::ofstream file(path);
    file << std::setw(4) << save << std::endl;
}

bool QuboSolution::isValid() const {
    return std::all_of(assignment_.begin(), assignment_.end(), [&](assignment_t u) { return u != kUNASSIGNED; });
}

QuboSolution::QuboSolution(const QUBO *qubo, const std::string &file) : qubo_(qubo) {
    int n = qubo_->getDim();
    assignment_ = std::vector<assignment_t>(n, kUNASSIGNED);

    std::ifstream f(file);
    nlohmann::json data = nlohmann::json::parse(f);

    for (int u : data[ASSIGN_0_NAME]) {
        assert(u < n);
        assert(assignment_[u] == kUNASSIGNED);
        assignment_[u] = 0;
    }

    for (int u : data[ASSIGN_1_NAME]) {
        assert(u < n);
        assert(assignment_[u] == kUNASSIGNED);
        assignment_[u] = 1;
    }
}

} // namespace sms