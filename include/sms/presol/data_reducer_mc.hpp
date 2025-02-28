#ifndef SMS_DATA_REDUCER_MC_HPP
#define SMS_DATA_REDUCER_MC_HPP

#include <cstdint>

#include "networkit/auxiliary/Timer.hpp"
#include "networkit/graph/Graph.hpp"
#include "tlx/container/d_ary_addressable_int_heap.hpp"

#include "sms/auxiliary/chrono_io.hpp"
#include "sms/instance/mc_solution.hpp"

namespace sms {

/**
 * @brief Flags for the reduction rules.
 * The higher the value the earlier the rule is checked.
 */
enum class RuleFlag : u_int16_t { // NOLINT(*-enum-size)
    kNONE = 0,
    kDegreeThree = 1 << 1, // might not delete edges
    kTriangle = 1 << 2,
    kSimilarNeighborhood = 1 << 3,
    kDominatingEdge = 1 << 4,
    kDegreeTwo = 1 << 5,
    // all rules below are weight stable, they do not introduce new edge weights
    kWeightStableClique = 1 << 6,    // nucleus is node with neighbors forming a clique
    kWeightStableThreePath = 1 << 7, // nucleus is degree two node with same weight neighbors
    kDegreeOne = 1 << 8,
    kDegreeZero = 1 << 9,
    kAll = (kDegreeZero << 1) - 1, // set all relevant bits to one
};

RuleFlag operator~(RuleFlag lhs);

RuleFlag operator|(RuleFlag lhs, RuleFlag rhs);

RuleFlag operator&(RuleFlag lhs, RuleFlag rhs);

RuleFlag highestPriorityRule(RuleFlag rf);

std::underlying_type_t<RuleFlag> to_underlying_type(RuleFlag rf);

struct NodeAggregation {
    NetworKit::node toAggregate;
    std::vector<NetworKit::node> others;
    bool toDifferent;
    bool negate1;
    bool negate2;
    // 0 deletion, 1 simple, 2 multi, 3 to smaller
    uint8_t type;

    // Vertex deletion
    explicit NodeAggregation(NetworKit::node toAggr) : toAggregate(toAggr), toDifferent(false), type(0) {
        assert(toAggr != NetworKit::none);
    }

    // Simple aggregation
    NodeAggregation(NetworKit::node toAggr, NetworKit::node n0, bool pToDiff)
        : toAggregate(toAggr), toDifferent(pToDiff), type(1) {
        others.push_back(n0);
        assert(n0 != NetworKit::none);
    }

    // Multi aggregation
    NodeAggregation(NetworKit::node toAggr, NetworKit::node n0, NetworKit::node n1, NetworKit::node n2, bool pToDiff,
                    bool neg1, bool neg2)
        : toAggregate(toAggr), toDifferent(pToDiff), negate1(neg1), negate2(neg2), type(2) {
        others.push_back(n0);
        others.push_back(n1);
        others.push_back(n2);
    }

    NodeAggregation(NetworKit::node toAggr, const std::vector<NetworKit::node> &relevantVertices)
        : toAggregate(toAggr), others(relevantVertices), type(3) {}

    /*
     * Given the partitioning of the neighbors, find partition of the aggregated node
     */
    bool reverse(const std::vector<bool> &partitions) const {
        if (partitions.size() != others.size()) {
            throw std::runtime_error("Partition size does not match number of neighbors!");
        }

        switch (type) {
            case 0:
                return toDifferent;
            case 1:
                return partitions[0] ^ toDifferent;
            case 2: {
                bool first = partitions[0] ^ (negate1 ^ partitions[1]);
                bool second = partitions[0] ^ (negate2 ^ partitions[2]);
                return (first & second) ^ toDifferent ^ partitions[0];
            }
            case 3: {
                auto ones = std::count(partitions.begin(), partitions.end(), true);
                return ones <= (std::ssize(partitions) / 2);
            }
            default:
                throw std::runtime_error("Unknown aggregation type!");
        }
    }
};

/*
 * From the stack of node aggregation and the solution to the transformed problem
 * Recover the original solution
 */
std::vector<Partition> recoverSolution(const std::vector<NodeAggregation> &aggr,
                                       const std::vector<Partition> &transformedSol);

NodeAggregation columnToNodeAggregation(NetworKit::node u, const std::array<NetworKit::node, 2> &separator,
                                        const std::array<uint8_t, 2> &column);

/*
 * Columns need to be color of node u when coloring of neighbors is in order
 * 0 0 0
 * 0 0 1
 * 0 1 0
 * 0 1 1
 */
NodeAggregation columnToNodeAggregation(NetworKit::node u, const std::array<NetworKit::node, 3> &neighbors,
                                        const std::array<uint8_t, 4> &column);

// lambda like object for the heap
template <class Type>
struct LargerInVector {
    explicit LargerInVector(const std::vector<Type> &vec) : vec_(vec) {}
    bool operator()(uint64_t x, uint64_t y) const noexcept { return (vec_)[x] > (vec_)[y]; }

private:
    const std::vector<Type> &vec_;
};

struct WeightedTriangle {
    NetworKit::node u;
    NetworKit::node v;
    NetworKit::node w;
    NetworKit::edgeweight uvWeight;
    NetworKit::edgeweight vwWeight;
    NetworKit::edgeweight wuWeight;
};

/**
 * @brief Data reducer for the MaxCut problem.
 */
class DataReducerMC {
public:
    /*
     * Creates a copy of the graph!
     */
    explicit DataReducerMC(const NetworKit::Graph &g, int emphasis = 1);

    // Needed because of the heap implementation
    DataReducerMC(DataReducerMC &&other) noexcept;

    // delete remaining default constructors
    DataReducerMC() = delete;
    DataReducerMC(const DataReducerMC &) = delete;

    DataReducerMC &operator=(const DataReducerMC &) = delete;
    DataReducerMC &operator=(DataReducerMC &&) = delete;

    NetworKit::edgeweight getOffset() const { return offset_; }

    const NetworKit::Graph &getGraph() const { return graph_; }

    const std::vector<NodeAggregation> &getHistory() const { return history_; }

    std::chrono::milliseconds getElapsedMilliSecs() const { return elapsed_; }

    // Automatically called by run!
    void initRun();

    /**
     * @brief Run the data reduction.
     * Fills the history stack with the applied rules.
     */
    void run();

    nlohmann::ordered_json getStats() const;

    std::vector<Partition> recoverSolution(const std::vector<Partition> &transformedSol) const;

    /**
     * Tries to perform the next aggregation step and returns the number of aggregations pushed to the history
     */
    int nextNodeAggregation();

    /**
     * @brief Apply the degree zero rule to a node.
     * @param node
     * @return 1 because exactly one vertex is deleted.
     */
    int degreeZeroApplication(NetworKit::node node);

    int degreeOneApplication(NetworKit::node node);

    int degreeTwoApplication(NetworKit::node node);

    int dominatingEdgeApplication(NetworKit::node u);

    int degreeThreeApplication(NetworKit::node node);

    int threePathApplication(NetworKit::node u);

    int cliqueApplication(NetworKit::node u);

    int gomoryHuApplication();

    int smallSeparatorApplication();

    void activateRule(NetworKit::node node, RuleFlag rule) {
        vertexStatus_[node] = (vertexStatus_[node] | rule);
        vertexHeap_.update(node);
    }

    void deactivateRule(NetworKit::node node, RuleFlag rule) {
        assert(std::popcount(to_underlying_type(rule)) == 1);
        vertexStatus_[node] = (vertexStatus_[node] & ~rule);
        vertexHeap_.update(node);
    }

    void deactivateAllRules(NetworKit::node node) {
        vertexStatus_[node] = RuleFlag::kNONE;
        vertexHeap_.update(node);
    }

    void smartActivateAll(NetworKit::node u);

    void removeVertex(NetworKit::node u) {
        graph_.removeNode(u);
        deactivateAllRules(u);
    }

    void addEdge(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        graph_.addEdge(u, v, weight);
        allWeightsAreOne_ &= (weight == 1.);
        weightOneNeighborhood_[u] = weightOneNeighborhood_[u] && (weight == 1.);
        weightOneNeighborhood_[v] = weightOneNeighborhood_[v] && (weight == 1.);
    }

    void contractNodes(NetworKit::node nodeToDelete, NetworKit::node nodeToKeep, bool toSamePartition);

    bool isRuleActive(NetworKit::node node, RuleFlag rule) const { return (vertexStatus_[node] & rule) == rule; }

    int triangleApplication(NetworKit::node u);

    int triangleApplicationCore(WeightedTriangle &wt);

    bool neighborhoodIsClique(NetworKit::node v, bool degreesChecked);

    /*
     * @brief Check if a rule is active and deactivate it.
     * Returns activity status before deactivation.
     */
    bool checkActiveAndDeactivate(NetworKit::node node, RuleFlag rule) {
        if (isRuleActive(node, rule)) {
            deactivateRule(node, rule);
            return true;
        }
        return false;
    }

    NetworKit::node getTop() {
        assert(!vertexHeap_.empty());
        return vertexHeap_.top();
    }

    NetworKit::node extractTop() { return vertexHeap_.extract_top(); }

    RuleFlag getDegreeFlag(NetworKit::node u) const;

    int removeSeparatedCliques();

    void turnOffNewPresolve() {
        emphasis_ = std::min(emphasis_, 1); // emphasis >= 2 uses vertex separators
        useDegreeThree_ = false;
        useTriangleThree_ = false;
        useTwins_ = false;
        useFullSeparatedCliques_ = false;
        useFastSolver_ = false;
    }

    void turnOffTriangles() {
        useTriangles_ = false;
        useTriangleThree_ = false;
    }

    void turnOffNewTriangles() { useTriangleThree_ = false; }

private:
    NetworKit::Graph graph_;
    std::vector<RuleFlag> vertexStatus_;
    tlx::DAryAddressableIntHeap<uint32_t, 2, LargerInVector<RuleFlag>> vertexHeap_;
    std::vector<bool> weightOneNeighborhood_; // all incident edges have weight 1
    bool allWeightsAreOne_;
    bool separateCliquesChecked_ = false;
    bool linearTimeSolveChecked_ = false;

    NetworKit::edgeweight offset_{0};
    std::vector<NetworKit::edgeweight> neighborWeights_;  // reusable array for the neighbor weights
    std::vector<NetworKit::edgeweight> neighborWeights2_; // reusable array for the neighbor weights
    std::vector<NetworKit::node> neighborMarks_;          // reusable array for the neighbor marks

    std::vector<NodeAggregation> history_;
    std::chrono::milliseconds elapsed_;

    // config
    int emphasis_;
    bool useDegreeThree_ = true;
    bool useTriangles_ = true;
    bool useTriangleThree_ = true;
    bool useTwins_ = true;
    bool useFullSeparatedCliques_ = true;
    bool useFastSolver_ = true;

    // statistics
    int posAlphaWithEdge_ = 0;
    int posAlphaWithoutEdge_ = 0;
    int negAlphaWithEdge_ = 0;
    int negAlphaWithoutEdge_ = 0;
    int twins_ = 0;
    int numCliques_ = 0;
    int numNearCliques_ = 0;
    int numCliquesVerticesRemoved_ = 0;
    int numSeparatedCliqueVerticesRemoved_ = 0;
    int numDominatingEdge_ = 0;
    int numRemovedFastSolver_ = 0;
    int numTriangles_ = 0;
    int numGomoryHuReductions_ = 0;
    int numSmallSeparatorReductions_ = 0;
    int degThreeReductions_ = 0;

    RuleFlag getNonDegreeFlags(NetworKit::node u) const;

    int oneSepApplication(const std::vector<NetworKit::node> &separatedSet, const std::vector<NetworKit::node> &sep);
    int twoSepApplication(const std::vector<NetworKit::node> &separatedSet, const std::vector<NetworKit::node> &sep);
    int threeSepApplication(const std::vector<NetworKit::node> &separatedSet, const std::vector<NetworKit::node> &sep);

    int similarNeighborhoodApplication(NetworKit::node u);
    int hashNeighborhoods();
    int cliqueReductionCore(NetworKit::node, unsigned int cliqueSize, unsigned int numExternalNodes);
    int nearCliqueReductionCore(NetworKit::node u, NetworKit::node partner, unsigned int cliqueSize,
                                unsigned int numExternalNodes);
};

} // namespace sms
#endif // SMS_DATA_REDUCER_MC_HPP
