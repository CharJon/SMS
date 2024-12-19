#ifndef SMS_BIPARTITE_HPP
#define SMS_BIPARTITE_HPP

#include "networkit/graph/Graph.hpp"

namespace sms {

class Bipartite {

public:
    explicit Bipartite(const NetworKit::Graph &g);

    void run(bool negativeWeightForcesSameSide = false);

    bool isBipartiteGraph() const {
        if (!hasRun_)
            throw std::runtime_error("Not yet run!");
        return isBipartite_;
    }

    const std::vector<uint8_t> &getPartition() const {
        if (!hasRun_)
            throw std::runtime_error("Not yet run!");
        return partition_;
    }

private:
    const NetworKit::Graph &graph_;
    std::vector<uint8_t> partition_;
    bool isBipartite_ = true;
    bool hasRun_ = false;
};

} // namespace sms

#endif // SMS_BIPARTITE_HPP
