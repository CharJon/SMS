#include "sms/graph/gomory_hu_tree.hpp"

#include "sms/graph/dinic.hpp"

namespace sms {

void GomoryHuTree::run() {

    for (auto e : graph_.edgeWeightRange()) {
        if (e.weight < 0)
            throw std::runtime_error("Only positive weights allowed!");
    }

    auto dinic = Dinic(graph_);

    for (NetworKit::node s : graph_.nodeRange()) {
        if (s == smallestNode_)
            continue;
        NetworKit::node t = representative_[s];
        dinic.run(s, t);
        auto flowValue = dinic.getFlowValue();
        auto sSide = dinic.getSourceSet();

        flow_[s] = flowValue;
        for (auto i : sSide) {
            if (i != s && representative_[i] == t)
                representative_[i] = s;
        }

        if (std::find(sSide.begin(), sSide.end(), representative_[t]) != sSide.end()) {
            representative_[s] = representative_[t];
            representative_[t] = s;
            flow_[s] = flow_[t];
            flow_[t] = flowValue;
        }
    }
    constructGomoryHuTree();
    run_ = true;
}

void GomoryHuTree::constructGomoryHuTree() {
    for (NetworKit::node s : graph_.nodeRange()) {
        if (s == smallestNode_)
            continue;
        auto x = representative_[s];
        gomoryHuTree_.addEdge(s, x, flow_[s]);
    }
}

const NetworKit::Graph &GomoryHuTree::getGomoryHuTree() const {
    if (!run_)
        throw std::runtime_error("The gomory hu tree hasn't been calculated yet.");
    return gomoryHuTree_;
}

// return minimum weight edge on s-t-path in gomory hu tree
NetworKit::WeightedEdge GomoryHuTree::minCutEdge(NetworKit::node s, NetworKit::node t) {
    s = contractedTo_[s];
    t = contractedTo_[t];
    if (!run_)
        throw std::runtime_error("The gomory hu tree hasn't been calculated yet.");
    std::vector<NetworKit::node> parent = std::vector<NetworKit::node>(graph_.upperNodeIdBound(), NetworKit::none);
    std::queue<NetworKit::node> queue;
    bool tReached = false;

    queue.push(s);

    while (!queue.empty()) {
        NetworKit::node current = queue.front();
        queue.pop();
        if (current == t)
            break;
        for (auto neighbor : gomoryHuTree_.neighborRange(current)) {
            if (parent[neighbor] != NetworKit::none)
                continue;
            queue.push(neighbor);
            parent[neighbor] = current;
            if (neighbor == t) {
                tReached = true;
                break;
            }
        }
        if (tReached)
            break;
    }

    if (parent[t] == NetworKit::none) {
        // Graph not connected
        return {NetworKit::none, NetworKit::none, 0};
    }

    double minCut = std::numeric_limits<double>::max();
    NetworKit::node current = t;
    NetworKit::WeightedEdge minEdge;

    while (current != s) {
        auto gomoryWeight = gomoryHuTree_.weight(current, parent[current]);
        minCut = std::min(minCut, gomoryWeight);
        if (minCut == gomoryWeight)
            minEdge = {current, parent[current], minCut};
        current = parent[current];
    }

    return minEdge;
}

void GomoryHuTree::contractNodesInGHT(NetworKit::node v, NetworKit::node w) {
    auto vInGHT = contractedTo_[v];
    auto wInGHT = contractedTo_[w];
    for (auto neighbor : gomoryHuTree_.neighborRange(vInGHT)) {
        if (neighbor != w) { // neighbor and w are not connected in the gomory hu tree
            gomoryHuTree_.addEdge(wInGHT, neighbor, gomoryHuTree_.weight(vInGHT, neighbor));
        }
    }
    for (auto node : contractions_[vInGHT]) {
        contractedTo_[node] = wInGHT;
        contractions_[wInGHT].push_back(node);
    }
    contractions_[vInGHT] = {};
    gomoryHuTree_.removeNode(vInGHT);
}

} // namespace sms