#include "sms/instance/maxcut.hpp"

#include "sms/auxiliary/math.hpp"
#include "sms/graph/graphs.hpp"

namespace sms {

MaxCut::MaxCut(const NetworKit::Graph &g, int shuffleVertices) : scalingFactor_(1.0), offset_(0.0) {
    auto res = compactGraph(g, shuffleVertices);
    graph_ = std::move(res.compactGraph);
    graph_.indexEdges();
    originalToNewNode_ = std::move(res.orig2compact);
    newToOriginalNode_ = std::move(res.compact2orig);

    for (auto e : graph_.edgeWeightRange()) {
        if (!isInteger(e.weight)) {
            integerSolutions_ = false;
            break;
        }
    }
}

double MaxCut::getSolutionValue(const std::vector<uint8_t> &solVector) const {
    assert(std::all_of(solVector.begin(), solVector.end(), [](auto i) { return i == 0 || i == 1; }));
    double res = 0;
    for (auto e : graph_.edgeWeightRange()) {
        res += (solVector[e.u] ^ solVector[e.v]) * e.weight;
    }
    return res + offset_;
}

nlohmann::ordered_json MaxCut::getInstanceInformation() {
    nlohmann::ordered_json j;

    j["num nodes"] = getNumberOfVertices();
    j["num edges"] = getNumberOfEdges();
    j["scaling factor"] = scalingFactor_;
    j["offset"] = offset_;

    return j;
}

NetworKit::edgeweight MaxCut::scale() {
    auto edgeWeightBasedDivisor = edgeWeightDivisor(graph_);
    if (edgeWeightBasedDivisor != 1.) {
        for (auto e : graph_.edgeWeightRange()) {
            auto newWeight = e.weight / edgeWeightBasedDivisor;
            graph_.setWeight(e.u, e.v, newWeight);
        }
        scalingFactor_ *= edgeWeightBasedDivisor;
    }

    auto degreeBasedDivisor = degreeBasedScaling(graph_);
    if (degreeBasedDivisor != 1.0) {
        for (auto e : graph_.edgeWeightRange()) {
            auto newWeight = e.weight / degreeBasedDivisor;
            graph_.setWeight(e.u, e.v, newWeight);
        }
        scalingFactor_ *= degreeBasedDivisor;
    }

    return edgeWeightBasedDivisor * degreeBasedDivisor;
}

void MaxCut::printInstanceInformation(std::ostream &out) {
    out << "---------- MaxCut instance statistics -----------------" << std::endl;
    auto j = getInstanceInformation();
    out << j.dump(4) << std::endl;
    out << "-------------------------------------------------------" << std::endl;
}

} // namespace sms