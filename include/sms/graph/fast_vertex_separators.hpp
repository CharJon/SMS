#ifndef SMS_FAST_VERTEX_SEPARATORS_HPP
#define SMS_FAST_VERTEX_SEPARATORS_HPP

#include "networkit/graph/Graph.hpp"

namespace sms {

class FastVertexSeparator {
public:
    explicit FastVertexSeparator(const NetworKit::Graph &g, unsigned int maxSeparatorSize = 3)
        : graph_(g),
          vertexInComponent_(graph_.upperNodeIdBound(), false),
          activeNodes_(graph_.upperNodeIdBound(), true),
          maxSeparatorSize_(maxSeparatorSize) {}

    void run(unsigned int maxSize);

    bool hasRun() const { return hasRun_; }

    size_t numberOfSeparators() const {
        assert(separators_.size() == separatedSet_.size());
        assertHasRun();
        return separators_.size();
    }

    const std::vector<std::vector<NetworKit::node>> &separators() const {
        assertHasRun();
        return separators_;
    }

    const std::vector<std::vector<NetworKit::node>> &separatedSets() const {
        assertHasRun();
        return separatedSet_;
    }

private:
    const NetworKit::Graph &graph_;
    bool hasRun_ = false;
    std::vector<std::vector<NetworKit::node>> separators_;
    std::vector<std::vector<NetworKit::node>> separatedSet_;
    std::vector<bool> vertexInComponent_;
    std::vector<bool> activeNodes_;
    unsigned int maxSeparatorSize_ = 3;
    unsigned int maxSize_ = 21;

    // find the smallest vertex separator, starting from node w
    void findSeparator(NetworKit::node);

    // finds the best vertex to add to the left side. Tries to strictly decrease size of separator
    size_t findBestNeighbor(std::vector<NetworKit::node>);

    int calcNewNeighborhoodSize(int oldNeighborhoodSize, NetworKit::node v);
    void updateNeighborhood(std::vector<NetworKit::node> &, NetworKit::node);
    void deactivateNodes(const std::vector<NetworKit::node> &);

    void assertHasRun() const {
        if (!hasRun_)
            throw std::runtime_error("Call run first!");
    }
};

} // namespace sms

#endif // SMS_FAST_VERTEX_SEPARATORS_HPP
