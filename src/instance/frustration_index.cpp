#include "sms/instance/frustration_index.hpp"

#include "sms/graph/graphs.hpp"

namespace sms {

FrustrationIndex::FrustrationIndex(const NetworKit::Graph &g) {
    assert(g.isWeighted());
    n_ = static_cast<int>(g.numberOfNodes());
    m_ = static_cast<int>(g.numberOfEdges());
    auto res = compactGraph(g);
    graph_ = std::move(res.compactGraph);
    originalToNewNode_ = std::move(res.orig2compact);
    newToOriginalNode_ = std::move(res.compact2orig);
}

double FrustrationIndex::getSolutionValue(const std::vector<uint8_t> &solVector) const {
    assert(std::all_of(solVector.begin(), solVector.end(), [](int i) { return i == 0 || i == 1; }));
    double res = 0;
    for (auto e : graph_.edgeWeightRange()) {
        if (e.weight == 1.0)
            res += (solVector[e.u] != solVector[e.v] ? 1.0 : 0.0);
        else if (e.weight == -1.0)
            res += (solVector[e.u] == solVector[e.v] ? 1.0 : 0.0);
    }

    return res;
}

bool FrustrationIndex::isConsistent() {
    for (auto e : graph_.edgeWeightRange()) {
        if (e.weight != 1.0 && e.weight != -1.0)
            return false;
    }
    return true;
}
nlohmann::ordered_json FrustrationIndex::getInstanceInformation() {
    nlohmann::ordered_json j;

    j["num nodes"] = n_;
    j["num edges"] = m_;

    int p = 0, n = 0;

    for (auto e : graph_.edgeWeightRange()) {
        if (e.weight == 1.0)
            p++;
        if (e.weight == -1.0)
            n++;
    }

    j["num negative edges"] = n;
    j["num positive edges"] = p;

    return j;
}
void FrustrationIndex::printInstanceInformation(std::ostream &out) {
    out << "---------- Frustration Index instance statistics -----------------" << std::endl;
    auto j = getInstanceInformation();
    out << j.dump(4) << std::endl;
    out << "------------------------------------------------------------------" << std::endl;
}

} // namespace sms
