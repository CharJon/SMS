#include "sms/graph/bipartite.hpp"

#include <cassert>

namespace sms {

constexpr uint8_t kUnassigned = 3;

Bipartite::Bipartite(const NetworKit::Graph &g) : graph_(g), partition_() {}

void Bipartite::run(bool negativeWeightForcesSameSide) {
    hasRun_ = true;
    isBipartite_ = true;
    partition_ = std::vector<uint8_t>(graph_.upperNodeIdBound(), kUnassigned);

    for (NetworKit::node v : graph_.nodeRange()) {
        if (partition_[v] != kUnassigned)
            continue;
        partition_[v] = 0;

        std::queue<NetworKit::node> queue;
        queue.push(v);

        while (not queue.empty()) {
            NetworKit::node w = queue.front();
            queue.pop();
            assert(partition_[w] != kUnassigned);

            if (!graph_.isWeighted()) {
                for (NetworKit::node x : graph_.neighborRange(w)) {
                    if (partition_[w] == partition_[x]) {
                        isBipartite_ = false;
                        return;
                    }

                    if (partition_[x] == kUnassigned) {
                        queue.push(x);
                        partition_[x] = 1 - partition_[w];
                    }
                }
            } else {
                for (auto [neighbor, weight] : graph_.weightNeighborRange(w)) {
                    if (weight < 0 && negativeWeightForcesSameSide) {
                        if (partition_[neighbor] == kUnassigned) {
                            queue.push(neighbor);
                            partition_[neighbor] = partition_[w];
                        } else if (partition_[w] != partition_[neighbor]) {
                            isBipartite_ = false;
                            return;
                        }
                    } else if (weight > 0 || !negativeWeightForcesSameSide) {
                        if (partition_[w] == partition_[neighbor]) {
                            isBipartite_ = false;
                            return;
                        }
                        if (partition_[neighbor] == kUnassigned) {
                            queue.push(neighbor);
                            partition_[neighbor] = 1 - partition_[w];
                        }
                    }
                }
            }
        }
    }
}

} // namespace sms
