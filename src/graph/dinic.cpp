#include "sms/graph/dinic.hpp"

namespace sms {

void Dinic::run(NetworKit::node source, NetworKit::node sink) {

    if (graph_.isEmpty()) {
        throw std::runtime_error("Graph is empty!");
    }
    if (!(graph_.hasNode(source) && graph_.hasNode(sink))) {
        throw std::runtime_error("Source or sink not in graph!");
    }
    if (graph_.isDirected()) {
        throw std::runtime_error("This Dinic implementation is for undirected graphs only!");
    }
    if (!graph_.isWeighted()) {
        throw std::runtime_error("Graph must contain weights!");
    }
    if (graph_.numberOfSelfLoops() > 0) {
        throw std::runtime_error("Graph contains Self Loops. Self Loops are not allowed!");
    }
    if (!graph_.hasEdgeIds()) {
        throw std::runtime_error("Graph needs to have edge IDs!");
    }

    if (source == sink) {
        maxFlowValue_ = std::numeric_limits<double>::max();
        hasRun_ = true;
        return;
    }
    source_ = source;
    sink_ = sink;
    maxFlowValue_ = 0.0;

    Dinic::resetFlow();

    bool done = false;
    while (!done) {
        done = !rankNodes() || !pathFinding();
    }
    hasRun_ = true;
}

bool Dinic::rankNodes() { // rank(v) = dist(v,t)

    nodeLevel_ = std::vector<int>(graph_.upperNodeIdBound(), -1);
    std::queue<NetworKit::node> q;

    nodeLevel_[sink_] = 0;
    q.push(sink_);

    while (!q.empty()) {
        NetworKit::node current = q.front();
        q.pop();
        if (current == source_)
            break;
        for (unsigned int i = 0; i < graph_.degree(current); ++i) {
            auto [neighbor, eid] = graph_.getIthNeighborWithId(current, i);
            auto weight = graph_.getIthNeighborWeight(current, i);
            if (nodeLevel_[neighbor] == -1 && residualCapacity(neighbor, current, weight, eid) > 0) {
                // only add those neighbors with rest capacity and which were not discovered yet
                nodeLevel_[neighbor] = nodeLevel_[current] + 1;
                q.push(neighbor);

                if (neighbor == source_) {
                    sourceLevel_ = nodeLevel_[neighbor];
                    return true;
                }
            }
        }
    }
    return false;
}

bool Dinic::pathFinding() {

    std::vector<bool> visited(graph_.upperNodeIdBound(), false);
    std::stack<NetworKit::node> stack;
    std::vector<NetworKit::node> path;

    stack.push(source_);
    visited[source_] = true;

    bool sinkReached = false;
    bool neighbourFound = false;

    while (!stack.empty()) {

        NetworKit::node current = stack.top();
        if (path.empty() || path.back() != current) {
            path.push_back(current);
        }

        // if s-t path was found, augment along it
        if (current == sink_) {

            double pathFlow = std::numeric_limits<double>::max();
            for (int i = 0; i < sourceLevel_; i++) {
                pathFlow =
                    std::min(pathFlow, residualCapacity(path[i], path[i + 1], graph_.weight(path[i], path[i + 1]),
                                                        graph_.edgeId(path[i], path[i + 1])));
            }
            for (int i = sourceLevel_; i > 0; i--) {
                NetworKit::edgeid eid = graph_.edgeId(path[i - 1], path[i]);
                NetworKit::edgeweight weight = graph_.weight(path[i - 1], path[i]);
                if (path[i - 1] < path[i]) {
                    flow_[eid] += pathFlow;
                    if (flow_[eid] == weight) {
                        current = path[i - 1];
                    }
                } else {
                    flow_[eid] -= pathFlow;
                    if (-flow_[eid] == weight) {
                        current = path[i - 1];
                    }
                }
            }
            maxFlowValue_ += pathFlow;

            sinkReached = true;
            NetworKit::node x = sink_;
            while (current != x) {
                stack.pop();
                path.pop_back();
                x = path.back();
            }
        }

        for (unsigned int i = 0; i < graph_.degree(current); ++i) {
            auto [neighbor, eid] = graph_.getIthNeighborWithId(current, i);
            auto weight = graph_.getIthNeighborWeight(current, i);
            if (residualCapacity(current, neighbor, weight, eid) > 0 && nodeLevel_[neighbor] == nodeLevel_[current] - 1
                && !visited[neighbor]) {
                stack.push(neighbor);
                neighbourFound = true;
                break;
            }
        }
        if (!neighbourFound) {
            visited[current] = true;
            stack.pop();
            path.pop_back();
        }
        neighbourFound = false;
    }
    return sinkReached;
}

std::vector<NetworKit::node> Dinic::getSourceSet() {
    if (!hasRun_)
        throw std::runtime_error("Error, run must be called first");

    // perform bfs from source
    std::vector<bool> visited(graph_.upperNodeIdBound(), false);
    std::vector<NetworKit::node> sourceSet;

    std::queue<NetworKit::node> queue;
    queue.push(source_);
    visited[source_] = true;
    while (!queue.empty()) {
        NetworKit::node u = queue.front();
        queue.pop();
        sourceSet.push_back(u);

        graph_.forNeighborsOf(
            u, [&](NetworKit::node, NetworKit::node v, NetworKit::edgeweight weight, NetworKit::edgeid eid) {
            if (!visited[v]
                && (residualCapacity(u, v, weight, eid) > 0 || residualCapacity(v, u, weight, eid) < weight)) {
                queue.push(v);
                visited[v] = true;
            }
        });
    }

    return sourceSet;
}

} // namespace sms
