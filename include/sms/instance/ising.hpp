#ifndef SMS_ISING_HPP
#define SMS_ISING_HPP

#include <optional>
#include <tuple>
#include <vector>

#include "nlohmann/json.hpp"

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"

namespace sms {

using NetworKit::Edge;
using NetworKit::Graph;
using NetworKit::node;

/*
 * Represents an Ising instance.
 * Objective is: Minimize \f$ -\sum_{i,j} q_{i,j} x_i x_j \f$
 * Where \f$ x_i \in \{-1, +1\} \f$.
 * Data is stored in a compact form using node id from 0 to n-1.
 */
class Ising {

public:
    explicit Ising(const Graph &g, double scalingFactor = 1, double offset = 0);

    virtual ~Ising() = default;

    int getNumberOfSpins() const { return n_; }

    int getNumberOfInteractions() const { return m_; }

    double getScalingFactor() const { return scalingFactor_; }

    double getOffset() const { return offset_; }

    double getSolutionValue(const std::vector<int> &solVector) const;

    const Graph &getGraph() const { return graph_; }

    node getOriginalNode(node newNode) const { return newToOriginalNode_.at(newNode); }

    node getNewNode(node originalNode) const { return originalToNewNode_.at(originalNode); }

    nlohmann::ordered_json getInstanceInformation();

    void printInstanceInformation(std::ostream &out);

protected:
    std::string name_;
    int n_;
    int m_;
    double scalingFactor_;
    double offset_;
    std::unordered_map<node, node> originalToNewNode_;
    std::vector<NetworKit::node> newToOriginalNode_;
    Graph graph_;
};

struct gridCoordinates {
    node layer;  // (optional) third dimension
    node row;    // second dimension
    node column; // first dimension

    bool operator==(const gridCoordinates &b) const {
        return (layer == b.layer) && (row == b.row) && (column == b.column);
    }
};

} // namespace sms

#endif // SMS_ISING_HPP
