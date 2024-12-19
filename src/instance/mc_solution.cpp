#include "sms/instance/mc_solution.hpp"

#include "nlohmann/json.hpp"

namespace sms {

nlohmann::json partitionToJson(const NetworKit::Graph &g, const std::vector<bool> &part) {
    nlohmann::json save;

    std::vector<NetworKit::node> part0nodes{};
    std::vector<NetworKit::node> part1nodes{};

    for (auto u : g.nodeRange()) {
        if (part[u]) {
            part0nodes.push_back(u);
        } else {
            part1nodes.push_back(u);
        }
    }

    save[kPartitionZeroName] = part0nodes;
    save[kPartitionOneName] = part1nodes;

    return save;
}

nlohmann::json partitionToJson(const std::vector<NetworKit::node> &existingNodes, const std::vector<bool> &part) {
    nlohmann::json save;

    std::vector<NetworKit::node> part0nodes{};
    std::vector<NetworKit::node> part1nodes{};

    for (auto u : existingNodes) {
        if (part[u]) {
            part0nodes.push_back(u);
        } else {
            part1nodes.push_back(u);
        }
    }

    save[kPartitionZeroName] = part0nodes;
    save[kPartitionOneName] = part1nodes;

    return save;
}

std::vector<bool> partitionFromJson(const nlohmann::json &j) {
    auto numNodes = j[kPartitionZeroName].size() + j[kPartitionOneName].size();
    std::vector<bool> part(numNodes, false);
    for (auto u : j[kPartitionOneName]) {
        part[u] = true;
    }
    return part;
}

NetworKit::edgeweight solutionValue(const NetworKit::Graph &g, const std::vector<Partition> &part) {
    NetworKit::edgeweight totalCutWeight = 0.0;
    for (auto e : g.edgeWeightRange()) {
        assert(part[e.u] != Partition::kUNASSIGNED);
        assert(part[e.v] != Partition::kUNASSIGNED);
        NetworKit::edgeweight value = e.weight * (part[e.u] ^ part[e.v]);
        totalCutWeight += value;
    }
    return totalCutWeight;
}

NetworKit::edgeweight solutionValue(const NetworKit::Graph &g, const std::vector<bool> &part) {
    NetworKit::edgeweight totalCutWeight = 0.0;
    for (auto e : g.edgeWeightRange()) {
        NetworKit::edgeweight value = e.weight * (part[e.u] ^ part[e.v]);
        totalCutWeight += value;
    }
    return totalCutWeight;
}

std::vector<NetworKit::node> allIds(const std::vector<bool> &values, bool match) {
    std::vector<NetworKit::node> ids;
    for (unsigned int i = 0; i < values.size(); ++i) {
        if (!(values[i] ^ match))
            ids.push_back(i);
    }
    return ids;
}

} // namespace sms