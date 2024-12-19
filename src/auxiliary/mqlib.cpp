#include "sms/auxiliary/mqlib.hpp"

std::vector<mqlib::Instance::InstanceTuple> mqlibEdgelist(const NetworKit::Graph &g) {
    std::vector<mqlib::Instance::InstanceTuple> result;
    result.reserve(g.numberOfEdges());
    for (auto e : g.edgeWeightRange()) {
        auto u = static_cast<int>(e.u) + 1; // move to one indexed
        auto v = static_cast<int>(e.v) + 1; // move to one indexed
        auto weight = e.weight;
        result.emplace_back(std::make_pair(u, v), weight);
    }
    return result;
}
