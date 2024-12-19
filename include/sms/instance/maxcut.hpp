#ifndef SMS_MAXCUT_HPP
#define SMS_MAXCUT_HPP

#include "nlohmann/json.hpp"

#include "networkit/graph/Graph.hpp"

namespace sms {

/**
 * Represents a MaxCut problem on a graph with n vertices.
 * Data is stored in a compact form using node id from 0 to n-1.
 */
class MaxCut {

public:
    /**
     * Creates a MaxCut instance from a graph. A compact copy of the graph is created and stored!
     * @param g The graph.
     * @param shuffleVertices Shall the vertices be shuffled, <= 0 means no shuffle
     */
    explicit MaxCut(const NetworKit::Graph &g, int shuffleVertices = 0);

    int getNumberOfVertices() const { return graph_.numberOfNodes(); }

    int getNumberOfEdges() const { return graph_.numberOfEdges(); }

    int maxOriginalNodeId() const { return *std::max_element(newToOriginalNode_.begin(), newToOriginalNode_.end()); }

    // a (not necessarily sorted) list of all original nodes
    const std::vector<NetworKit::node> &originalNodes() const { return newToOriginalNode_; }

    void setScalingFactor(double scaling) { scalingFactor_ = scaling; }

    double getScalingFactor() const { return scalingFactor_; }

    /*
     * Scales the edge weights as small as possible, but keeps solution values integer (if they were before scaling).
     */
    NetworKit::edgeweight scale();

    void setOffset(double offset) { offset_ = offset; }

    double getOffset() const { return offset_; }

    double getSolutionValue(const std::vector<uint8_t> &solVector) const;

    const NetworKit::Graph &getGraph() const { return graph_; }

    NetworKit::node getOriginalNode(NetworKit::node newNode) const { return newToOriginalNode_.at(newNode); }

    NetworKit::node getNewNode(NetworKit::node originalNode) const { return originalToNewNode_.at(originalNode); }

    bool isOriginalNode(NetworKit::node node) const {
        return originalToNewNode_.find(node) != originalToNewNode_.end();
    }

    const std::string &getName() const { return name_; }

    nlohmann::ordered_json getInstanceInformation();

    void printInstanceInformation(std::ostream &out);

    bool allSolutionsHaveIntegerWeight() const { return integerSolutions_; }

private:
    std::string name_;
    double scalingFactor_;
    double offset_;
    NetworKit::Graph graph_;
    std::unordered_map<NetworKit::node, NetworKit::node> originalToNewNode_;
    std::vector<NetworKit::node> newToOriginalNode_;
    bool integerSolutions_ = true;
};

} // namespace sms

#endif // SMS_MAXCUT_HPP
