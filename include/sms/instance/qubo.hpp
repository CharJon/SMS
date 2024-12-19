#ifndef SMS_QUBO_HPP
#define SMS_QUBO_HPP

#include <string>
#include <vector>

#include "nlohmann/json.hpp"

#include "sms/io/io.hpp"

namespace sms {
/*
 * Represents a QUBO instance with dimension d. QUBOs are minimization problems here.
 * Data is stored in a compact form using variable ids from 0 to d-1.
 */
class QUBO {

public:
    explicit QUBO(int dim);

    QUBO(const double *q, int dim, double offset = 0.0, double scalingFactor = 1.0);

    void setValue(unsigned int i, unsigned int j, double value);

    double getValue(unsigned int i, unsigned int j) const;

    /***
     * @brief Get solution value for core data (scaling and offset are ignored)
     */
    template <typename T>
    double getSolutionValue(const T *solVector) const {
        double res = 0;
        for (int i = 0; i < dim_; i++) {
            for (int j = 0; j < dim_; j++) {
                res += matrix_[i * dim_ + j] * solVector[i] * solVector[j];
            }
        }
        return res;
    }

    /***
     * @brief Get the original solution value, scaling and offset are applied
     */
    template <typename T>
    double getOriginalSolutionValue(const T *solVector) const {
        return getSolutionValue(solVector) * scalingFactor_ + offset_;
    }

    double getSolutionValueBLAS(const double *solVector) const;

    double getOriginalSolutionValueBLAS(const double *solVector) const;

    int getDim() const;

    bool isValid() const;

    void resetToZero();

    void fillAll(double value);

    double getScalingFactor() const { return scalingFactor_; }

    double getOffset() const { return offset_; }

    void setScalingFactor(double scalingFactor) { scalingFactor_ = scalingFactor; }

    void setOffset(double offset) { offset_ = offset; }

    nlohmann::ordered_json getInstanceInformation();

    void printInstanceInformation(std::ostream &out);

private:
    std::string name_;
    int dim_;
    std::vector<double> matrix_;

    double scalingFactor_;
    double offset_;

    unsigned posInMatrix(unsigned int i, unsigned int j) const;
};

QUBO bqFileToQubo(const std::string &path);

} // namespace sms

#endif // SMS_QUBO_HPP
