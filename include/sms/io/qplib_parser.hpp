#ifndef SMS_QPLIB_PARSER_HPP
#define SMS_QPLIB_PARSER_HPP

#include <string>
#include <vector>

#include "sms/instance/qubo.hpp"
#include "sms/io/io.hpp"

namespace sms {
/***
 * See wiki/qplib_file_format.md for the file format specification.
 * For simplicity we add all linear terms to the diagonal.
 * Note: The quadratic part has a 0.5 scaling per qpblib specification!
 */
class QPLibParser {
public:
    explicit QPLibParser(const std::string &path) : instance_(0) {
        unparsedContent_ = fileToVector(path);
        removeComments();
        parseHeader(); // sets dim_
        instance_ = QUBO(dim_);
        parseQuadraticObjContent();
        parseLinearObjHeader();
        parseLinearObjContent();
        parseRemainder();
        instance_.setOffset(constantOffset_);
        instance_.setScalingFactor(isMaximize_ ? -1 : 1); // quadratic part is 1/2 x^T Q x
    }

    /***
     * @return Copy of internally build instance.
     */
    QUBO getInstance();

private:
    std::vector<std::string> unparsedContent_;

    QUBO instance_;
    int numberNoDefaultDiagonal_;
    double linearObjCoefDefault_;
    bool isMaximize_;
    int dim_;
    int numberOfTerms_;
    double constantOffset_;

    void removeComments();

    void parseHeader();

    void parseQuadraticObjContent();

    void parseLinearObjHeader();

    void parseLinearObjContent();

    void parseRemainder();
};

} // namespace sms

#endif // SMS_QPLIB_PARSER_HPP
