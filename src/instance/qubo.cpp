#include "sms/instance/qubo.hpp"

#include <algorithm>
#include <cassert>

#include "sms/auxiliary/math.hpp"
#include "sms/io/io.hpp"

namespace sms {

unsigned int QUBO::posInMatrix(unsigned int i, unsigned int j) const {
    return (i * dim_ + j);
}

QUBO::QUBO(int dim) : dim_(dim) {
    dim_ = dim;
    auto sizeDim = static_cast<size_t>(dim);
    matrix_ = std::vector<double>(sizeDim * sizeDim, 0.0);
    resetToZero();

    scalingFactor_ = 1;
    offset_ = 0;
}

QUBO::QUBO(const double *q, int dim, double offset, double scalingFactor) : dim_(dim) {
    auto sizeDim = static_cast<size_t>(dim);
    matrix_ = std::vector<double>(sizeDim * sizeDim, 0.0);

    for (int i = 0; i < dim * dim; i++) {
        matrix_[i] = q[i];
    }

    scalingFactor_ = scalingFactor;
    offset_ = offset;
}

void QUBO::setValue(unsigned int i, unsigned int j, double value) {
    assert(i < static_cast<unsigned int>(dim_) && j < static_cast<unsigned int>(dim_));
    matrix_[posInMatrix(i, j)] = value;
}

double QUBO::getValue(unsigned int i, unsigned int j) const {
    assert(i < static_cast<unsigned int>(dim_) && j < static_cast<unsigned int>(dim_));
    return matrix_[posInMatrix(i, j)];
}

double QUBO::getSolutionValueBLAS(const double *solVector) const {
    return quboEvaluation(&matrix_[0], solVector, dim_, 1.0);
}

double QUBO::getOriginalSolutionValueBLAS(const double *solVector) const {
    return quboEvaluation(&matrix_[0], solVector, dim_, 1.0) * scalingFactor_ + offset_;
}

int QUBO::getDim() const {
    return dim_;
}

bool QUBO::isValid() const {
    return (dim_ * dim_) == static_cast<int>(matrix_.size());
}

void QUBO::resetToZero() {
    fillAll(0);
}

void QUBO::fillAll(double value) {
    std::fill_n(matrix_.begin(), dim_ * dim_, value);
}

nlohmann::ordered_json QUBO::getInstanceInformation() {
    nlohmann::ordered_json j;

    j["dimension"] = dim_;
    j["scaling factor"] = scalingFactor_;
    j["offset"] = offset_;

    return j;
}

void QUBO::printInstanceInformation(std::ostream &out) {
    out << "---------- MaxCut instance statistics -----------------" << std::endl;
    auto j = getInstanceInformation();
    out << j.dump(4) << std::endl;
    out << "-------------------------------------------------------" << std::endl;
}

QUBO bqFileToQubo(const std::string &path) {
    auto bq = fileToBqObj(path);

    QUBO qubo(static_cast<int>(bq.nodes));

    for (auto [a, b, w] : bq.edge_list) {
        qubo.setValue(a - 1, b - 1, w);
    }

    return qubo;
}

} // namespace sms