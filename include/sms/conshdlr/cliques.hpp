#ifndef SMS_CLIQUES_HPP
#define SMS_CLIQUES_HPP

#include <cmath>
#include <optional>

#include "networkit/graph/Graph.hpp"

namespace sms {

class CliqueSeparator {

public:
    explicit CliqueSeparator(const NetworKit::Graph &originalGraph);

    void updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w);

    void addClique(const std::vector<NetworKit::node> &clique) { maximalCliques_.push_back(clique); }

    const std::vector<NetworKit::node> &getClique(unsigned int pos) const;

    size_t numCliques() const { return maximalCliques_.size(); }

    size_t numLargeCliques() const {
        return std::count_if(maximalCliques_.begin(), maximalCliques_.end(),
                             [this](auto c) { return std::size(c) > smallCliqueSize_; });
    }

    /*
     * Main function to check for violated constraints
     * Checks one clique and cliques contained inside this clique for violated constraints
     * Split is empty, if no violated constraint is found
     */
    std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>>
    checkViolation(const std::vector<NetworKit::node> &clique, int lowerSubCliqueSize = 5,
                   int upperSubCliqueSize = std::numeric_limits<int>::max(), bool recursive = true) const;

    // Checks one clique for a violated constraint
    // Split is empty, if no violated constraint is found
    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
    checkSimpleViolation(std::vector<NetworKit::node> clique) const;

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
    checkBigClique(const std::vector<NetworKit::node> &clique) const;

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
    checkBigOddClique(const std::vector<NetworKit::node> &clique) const;

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
    checkBigOddCliqueSimple(const std::vector<NetworKit::node> &clique) const;

    void sortCliquesBySize();

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
    checkTriPartition(const std::vector<NetworKit::node> &) const;

    std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>>
    checkSubCliques(const std::vector<NetworKit::node> &clique, int subCliqueSize, int lowerSubCliqueSize = 5,
                    int upperSubCliqueSize = std::numeric_limits<int>::max(), bool recursive = true) const;

private:
    NetworKit::Graph lpWeightedGraph_;
    std::vector<std::vector<NetworKit::node>> maximalCliques_;

    // config
    unsigned int smallCliqueSize_ = 24;
    unsigned int enumerationLimit_ = 12;
    bool fromSmallToBig_ = true;
};

} // namespace sms

#endif // SMS_CLIQUES_HPP
