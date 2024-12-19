#include "sms/io/qplib_parser.hpp"

namespace sms {

void QPLibParser::removeComments() {
    std::vector<std::string> res;
    for (auto line : unparsedContent_) {
        if ((line[0] != '!') && (line[0] != '#') && (line[0] != '%'))
            res.push_back(line);
    }

    unparsedContent_ = res;
}

void QPLibParser::parseHeader() {
    std::string type = unparsedContent_[1];
    std::string sense = unparsedContent_[2];

    auto t = splitString(type, ' ');

    if ((t[0] != "QBB") && (t[0] != "QBN")) {
        throw std::runtime_error("Type " + t[0] + " not supported!");
    }
    auto s = splitString(sense, ' ');

    if (s[0] == "maximize") {
        isMaximize_ = true;
    } else if (s[0] == "minimize") {
        isMaximize_ = false;
    } else {
        throw std::runtime_error("Objective type " + s[0] + " not supported!");
    }

    std::string dimLine = unparsedContent_[3];
    std::string termsLine = unparsedContent_[4];

    dim_ = strictStringToInt(splitString(dimLine, ' ')[0]);
    numberOfTerms_ = strictStringToInt(splitString(termsLine, ' ')[0]);

    unparsedContent_.erase(unparsedContent_.begin(), unparsedContent_.begin() + 5);
}

void QPLibParser::parseQuadraticObjContent() {

    for (int i = 0; i < numberOfTerms_; i++) {
        auto line = unparsedContent_[i];

        auto split = splitString(line, ' ');

        if (split.size() < 3)
            throw std::runtime_error("Matrix line has wrong format: " + line);

        // QPLib uses 1-based indexing
        int rowIdx = strictStringToInt(split[0]) - 1;
        int colIdx = strictStringToInt(split[1]) - 1;
        double value = strictStringToDouble(split[2]);

        if (colIdx > rowIdx)
            throw std::runtime_error("Entry (" + line + ") is not in lower triangle");

        instance_.setValue(rowIdx, colIdx, isMaximize_ ? -0.5 * value : 0.5 * value);
    }

    unparsedContent_.erase(unparsedContent_.begin(), unparsedContent_.begin() + numberOfTerms_);
}

void QPLibParser::parseLinearObjHeader() {

    auto defaultLine = unparsedContent_[0];
    auto nonDefaultLine = unparsedContent_[1];

    linearObjCoefDefault_ = strictStringToDouble(splitString(defaultLine, ' ')[0]);
    numberNoDefaultDiagonal_ = strictStringToInt(splitString(nonDefaultLine, ' ')[0]);

    unparsedContent_.erase(unparsedContent_.begin(), unparsedContent_.begin() + 2);
}

void QPLibParser::parseLinearObjContent() {
    // fill diagonal with default value times 2 because of 0.5 scaling
    auto scaledDefault = linearObjCoefDefault_ * (isMaximize_ ? -1 : 1);
    for (int i = 0; i < dim_; i++) {
        assert(instance_.getValue(i, i) == 0.0);
        auto newValue = instance_.getValue(i, i) + scaledDefault;
        instance_.setValue(i, i, newValue);
    }

    for (int i = 0; i < numberNoDefaultDiagonal_; i++) {
        auto line = unparsedContent_[i];
        auto split = splitString(line, ' ');
        if (split.size() < 2)
            throw std::runtime_error("Diagonal line has wrong format: " + line);

        auto rowIdx = strictStringToInt(split[0]) - 1; // from 1 indexed to zero indexed
        auto rowValue = strictStringToDouble(split[1]);
        auto scaledRowValue = rowValue * (isMaximize_ ? -1 : 1);
        auto currentValue = instance_.getValue(rowIdx, rowIdx);
        auto newValue = currentValue - scaledDefault + scaledRowValue;
        instance_.setValue(rowIdx, rowIdx, newValue);
    }

    unparsedContent_.erase(unparsedContent_.begin(), unparsedContent_.begin() + numberNoDefaultDiagonal_);
}

QUBO QPLibParser::getInstance() {
    return instance_;
}

void QPLibParser::parseRemainder() {
    auto offsetLine = unparsedContent_[0];

    constantOffset_ = strictStringToDouble(splitString(offsetLine, ' ')[0]);
    unparsedContent_.erase(unparsedContent_.begin(), unparsedContent_.begin() + 1);
}

} // namespace sms
