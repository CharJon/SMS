#include "sms/presol/data_reducer_mc.hpp"

#include <cstdlib>

#include "sms/auxiliary/chrono_io.hpp"
#include "sms/auxiliary/math.hpp"
#include "sms/graph/fast_vertex_separators.hpp"
#include "sms/graph/gomory_hu_tree.hpp"
#include "sms/graph/graphs.hpp"
#include "sms/solver/enumeration_solver.hpp"
#include "sms/solver/special_structures_solver.hpp"

namespace sms {

RuleFlag operator~(const RuleFlag operand) {
    return static_cast<RuleFlag>(~to_underlying_type(operand));
}

RuleFlag operator|(const RuleFlag lhs, const RuleFlag rhs) {
    return static_cast<RuleFlag>(to_underlying_type(lhs) | to_underlying_type(rhs));
}

RuleFlag operator&(const RuleFlag lhs, const RuleFlag rhs) {
    return static_cast<RuleFlag>(to_underlying_type(lhs) & to_underlying_type(rhs));
}

/*
 * Extracts the rule with the highest priority from a given RuleFlag.
 */
RuleFlag highestPriorityRule(const RuleFlag rf) {
    const auto leadZero = std::countl_zero(to_underlying_type(rf));
    const auto numShifts = std::max(15 - leadZero, 0);
    return static_cast<RuleFlag>(static_cast<std::underlying_type_t<RuleFlag>>(1) << numShifts);
}

std::underlying_type_t<RuleFlag> to_underlying_type(RuleFlag rf) {
    return static_cast<std::underlying_type_t<RuleFlag>>(rf);
}

NetworKit::edgeweight absOutWeight(const NetworKit::Graph &graph, NetworKit::node u) {
    NetworKit::edgeweight sum = 0;

    for (auto [neighbor, weight] : graph.weightNeighborRange(u)) {
        sum += std::abs(weight);
    }

    return sum;
}

int cutSizePerfectSplit(unsigned int numNodes) {
    int numSmallerSide = static_cast<int>(numNodes / 2);             // div 2 floor
    int numLargerSide = static_cast<int>(numNodes) - numSmallerSide; // equal or larger by one
    return numSmallerSide * numLargerSide;
}

// for a given triangle, check if not cutting the first edge is always better than cutting it
bool canContractThreeNegative(NetworKit::edgeweight w01, NetworKit::edgeweight w02, NetworKit::edgeweight w12,
                              NetworKit::edgeweight absOutV0, NetworKit::edgeweight absOutV1,
                              NetworKit::edgeweight absOutV2) {
    assert((w01 < 0) && (w02 < 0) && (w12 < 0));

    // e_{0,1} and e_{0,2} are cut
    auto lhs1 = -w01 - w02;
    bool conditionCut1 = lhs1 >= absOutV0 - std::abs(w01) - std::abs(w02);
    bool conditionCut23 = lhs1 >= absOutV1 - std::abs(w01) - std::abs(w12) + absOutV2 - std::abs(w02) - std::abs(w12);
    // e_{0,1} and e_{1,2} are cut
    auto lhs2 = -w01 - w12;
    bool conditionCut2 = lhs2 >= absOutV1 - std::abs(w01) - std::abs(w12);
    bool conditionCut13 = lhs2 >= absOutV0 - std::abs(w01) - std::abs(w02) + absOutV2 - std::abs(w02) - std::abs(w12);

    return (conditionCut1 || conditionCut23) && (conditionCut2 || conditionCut13);
}

DataReducerMC::DataReducerMC(const NetworKit::Graph &g, int emphasis)
    : graph_(g),
      vertexStatus_(g.upperNodeIdBound(), RuleFlag::kAll),
      vertexHeap_(LargerInVector<sms::RuleFlag>{vertexStatus_}),
      weightOneNeighborhood_(g.upperNodeIdBound(), true),
      neighborWeights_(g.upperNodeIdBound(), 0.0),
      neighborWeights2_(g.upperNodeIdBound(), 0.0),
      neighborMarks_(g.upperNodeIdBound(), NetworKit::none),
      emphasis_(emphasis) {
    if (graph_.isDirected())
        throw std::invalid_argument("DataReducerMC only supports undirected graphs.");
    if (!graph_.isWeighted())
        throw std::invalid_argument("DataReducerMC only supports weighted graphs.");
    if (graph_.upperNodeIdBound() != graph_.numberOfNodes())
        throw std::invalid_argument("DataReducerMC only supports compact graphs.");
    if (graph_.numberOfSelfLoops() > 0)
        throw std::invalid_argument("DataReducerMC does not support self loops.");
}

void DataReducerMC::initRun() {
    vertexHeap_.clear();
    allWeightsAreOne_ = true;
    for (auto e : graph_.edgeWeightRange()) {
        neighborWeights_[e.u] += e.weight < 0; // count number of incident negative edges
        neighborWeights_[e.v] += e.weight < 0; // count number of incident negative edges
        weightOneNeighborhood_[e.u] = weightOneNeighborhood_[e.u] && (e.weight == 1.0);
        weightOneNeighborhood_[e.v] = weightOneNeighborhood_[e.v] && (e.weight == 1.0);
        allWeightsAreOne_ &= (e.weight == 1.0);
    }
    std::fill(neighborWeights_.begin(), neighborWeights_.end(), 0.0); // reset vector
    separateCliquesChecked_ = !allWeightsAreOne_;
    linearTimeSolveChecked_ = false;

    for (auto node : graph_.nodeRange()) {
        RuleFlag rules = getDegreeFlag(node);
        rules = static_cast<RuleFlag>(rules | getNonDegreeFlags(node));
        vertexStatus_[node] = rules;
        if (graph_.degree(node) == 0) {
            assert(vertexStatus_[node] & RuleFlag::kDegreeZero);
        }
        vertexHeap_.push(node);
    }
    assert(vertexHeap_.sanity_check());
}

DataReducerMC::DataReducerMC(DataReducerMC &&other) noexcept
    : graph_(std::move(other.graph_)),
      vertexStatus_(std::move(other.vertexStatus_)),
      vertexHeap_(LargerInVector<sms::RuleFlag>{vertexStatus_}),
      weightOneNeighborhood_(std::move(other.weightOneNeighborhood_)),
      allWeightsAreOne_(other.allWeightsAreOne_),
      offset_(other.offset_),
      neighborWeights_(std::move(other.neighborWeights_)),
      neighborWeights2_(std::move(other.neighborWeights2_)),
      neighborMarks_(std::move(other.neighborMarks_)),
      history_(std::move(other.history_)),
      elapsed_(other.elapsed_),
      emphasis_(other.emphasis_),
      useDegreeThree_(other.useDegreeThree_),
      useTriangles_(other.useTriangles_),
      useTriangleThree_(other.useTriangleThree_),
      useTwins_(other.useTwins_),
      useFullSeparatedCliques_(other.useFullSeparatedCliques_),
      useFastSolver_(other.useFastSolver_) {}

void DataReducerMC::run() {
    Aux::StartedTimer t;
    initRun();

    // remove edges with a weight of zero
    std::vector<NetworKit::Edge> edgesToRemove;
    for (auto e : graph_.edgeWeightRange()) {
        if (e.weight == 0.0) {
            edgesToRemove.push_back(e);
        }
    }
    for (auto e : edgesToRemove) {
        graph_.removeEdge(e.u, e.v);

        if (graph_.degree(e.u) == 0) {
            vertexStatus_[e.u] = RuleFlag::kDegreeZero;
            vertexHeap_.update(e.u);
        }

        if (graph_.degree(e.v) == 0) {
            vertexStatus_[e.u] = RuleFlag::kDegreeZero;
            vertexHeap_.update(e.u);
        }
    }

    int numReductions = emphasis_; // if emphasis <= 0, no reductions are performed
    while (numReductions > 0) {
        while (numReductions > 0) {
            while (numReductions > 0) {
                numReductions = emphasis_ > 0 ? nextNodeAggregation() : 0;
            }
            numReductions = emphasis_ > 0 ? hashNeighborhoods() : 0;
        }
        numReductions = emphasis_ >= 2 ? smallSeparatorApplication() : 0;
        numReductions += emphasis_ >= 4 ? gomoryHuApplication() : 0;
    }

    t.stop();
    elapsed_ = std::chrono::duration_cast<std::chrono::milliseconds>(t.stopTime() - t.startTime());
}

nlohmann::ordered_json DataReducerMC::getStats() const {
    nlohmann::ordered_json j;
    j["emphasis"] = emphasis_;
    j["#cliques"] = numCliques_;
    j["#near cliques"] = numNearCliques_;
    j["#cliques vertices removed"] = numCliquesVerticesRemoved_;
    j["#separated clique vertices removed"] = numSeparatedCliqueVerticesRemoved_;
    j["#dominating edges"] = numDominatingEdge_;
    j["#triangles"] = numTriangles_;
    j["#same neighborhood"] = posAlphaWithEdge_ + negAlphaWithEdge_ + posAlphaWithoutEdge_ + negAlphaWithoutEdge_;
    j["#twins"] = twins_;
    j["#gomory hu"] = numGomoryHuReductions_;
    j["#small separator"] = numSmallSeparatorReductions_;
    j["#removed by fast solver"] = numRemovedFastSolver_;
    j["#degree three"] = degThreeReductions_;
    j["elapsed [ms]"] = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_).count();
    return j;
}

int DataReducerMC::nextNodeAggregation() {
    auto curNode = getTop();
    while (vertexStatus_[curNode] != RuleFlag::kNONE) {
        RuleFlag highestActive = highestPriorityRule(vertexStatus_[curNode]);
#ifndef NDEBUG
        auto nodeRange = graph_.nodeRange();
        for (auto u : nodeRange) {
            assert(highestActive == RuleFlag::kDegreeZero || graph_.degree(u) > 0);
        }
#endif
        assert(highestActive == RuleFlag::kDegreeZero
               || std::all_of(nodeRange.begin(), nodeRange.end(), [&](auto v) { return graph_.degree(v) > 0; }));
        int numRemovedVertices = 0;

        // First check if we enter non weight stable phase
        if (!separateCliquesChecked_ && (highestActive < RuleFlag::kWeightStableClique)) {
            separateCliquesChecked_ = true;
            return removeSeparatedCliques();
        }
        if (!linearTimeSolveChecked_ && (highestActive < RuleFlag::kWeightStableClique) && useFastSolver_) {
            linearTimeSolveChecked_ = true;
            auto compGraph = compactGraph(graph_);
            SpecialStructuresSolver solver(compGraph.compactGraph);
            if (solver.applicable()) {
                solver.run();
                auto sol = solver.getBestSolution();
                // remove one by one and map to arbitrary root 0
                auto origRoot = compGraph.compact2orig[0];
                for (unsigned int u = 1; u < compGraph.compactGraph.numberOfNodes(); u++) {
                    auto origU = compGraph.compact2orig[u];
                    contractNodes(origU, origRoot, sol[u] == sol[0]);
                }
                numRemovedFastSolver_ += compGraph.compactGraph.numberOfEdges();
                degreeZeroApplication(origRoot);
                assert(graph_.numberOfNodes() == 0);
                return compGraph.compactGraph.numberOfNodes();
            }
        }

        // In general
        switch (highestActive) {
            case RuleFlag::kDegreeZero:
                deactivateRule(curNode, RuleFlag::kDegreeZero);
                numRemovedVertices = degreeZeroApplication(curNode);
                break;
            case RuleFlag::kDegreeOne:
                deactivateRule(curNode, RuleFlag::kDegreeOne);
                numRemovedVertices = degreeOneApplication(curNode);
                break;
            case RuleFlag::kWeightStableThreePath:
                deactivateRule(curNode, RuleFlag::kWeightStableThreePath);
                numRemovedVertices = threePathApplication(curNode);
                break;
            case RuleFlag::kDegreeTwo:
                deactivateRule(curNode, RuleFlag::kDegreeTwo);
                numRemovedVertices = degreeTwoApplication(curNode);
                break;
            case RuleFlag::kDominatingEdge:
                deactivateRule(curNode, RuleFlag::kDominatingEdge);
                numRemovedVertices = dominatingEdgeApplication(curNode);
                break;
            case RuleFlag::kSimilarNeighborhood:
                deactivateRule(curNode, RuleFlag::kSimilarNeighborhood);
                numRemovedVertices = similarNeighborhoodApplication(curNode);
                break;
            case RuleFlag::kDegreeThree:
                deactivateRule(curNode, RuleFlag::kDegreeThree);
                numRemovedVertices = degreeThreeApplication(curNode);
                break;
            case RuleFlag::kTriangle:
                deactivateRule(curNode, RuleFlag::kTriangle);
                numRemovedVertices = triangleApplication(curNode);
                break;
            case RuleFlag::kWeightStableClique:
                deactivateRule(curNode, RuleFlag::kWeightStableClique);
                numRemovedVertices = cliqueApplication(curNode);
                numCliquesVerticesRemoved_ += numRemovedVertices;
                break;
            default:
                deactivateAllRules(curNode);
        }
        if (numRemovedVertices > 0)
            return numRemovedVertices;
        curNode = getTop();
    }

    return 0;
}

int DataReducerMC::hashNeighborhoods() {
    std::vector<NetworKit::node> hashes(graph_.upperNodeIdBound(), 0);

    for (auto u : graph_.nodeRange()) {
        NetworKit::node hash = 0;
        for (auto neighbor : graph_.neighborRange(u)) {
            hash ^= neighbor;
        }
        hashes[u] = hash;
    }

    auto hashFunc = [hashes](const NetworKit::node u) { return hashes[u]; };
    auto hashBuckets = std::unordered_set<NetworKit::node, decltype(hashFunc)>(2 * graph_.numberOfNodes(), hashFunc);

    for (auto u : graph_.nodeRange()) {
        hashBuckets.insert(u);
    }

    size_t max = 0;
    auto cnt = 0;
    for (size_t i = 0; i < hashBuckets.bucket_count(); i++) {
        auto curBucketSize = hashBuckets.bucket_size(i);
        max = std::max(max, curBucketSize);
        // for all pairs of nodes in bucket, do sth
        if (curBucketSize > 1) {
            for (auto firstIt = hashBuckets.begin(i); firstIt != hashBuckets.end(i); ++firstIt) {
                NetworKit::node u = *firstIt;
                if (!graph_.hasNode(u))
                    continue;
                auto secondIt = firstIt;
                ++secondIt;
                for (; secondIt != hashBuckets.end(i); ++secondIt) {
                    NetworKit::node v = *secondIt;
                    if (!graph_.hasNode(v))
                        continue;
                    auto alpha = neighborhoodAlpha(graph_, u, v, false, neighborWeights_);
                    if (alpha > 0.0) {
                        posAlphaWithoutEdge_ += 1;
                        contractNodes(v, u, true);
                        cnt += 1;
                    } else if (alpha < 0.0) {
                        negAlphaWithoutEdge_ += 1;
                        contractNodes(v, u, false);
                        cnt += 1;
                    }
                }
            }
        }
    }

    return cnt;
}

int DataReducerMC::degreeZeroApplication(NetworKit::node node) {
    if (graph_.degree(node) != 0)
        return 0;
    removeVertex(node);
    history_.emplace_back(node);
    return 1;
}

int DataReducerMC::degreeOneApplication(NetworKit::node node) {
    if (graph_.degree(node) != 1)
        return 0;
    auto [neighbor, neighborWeight] = graph_.getIthNeighborWithWeight(node, 0);
    contractNodes(node, neighbor, neighborWeight <= 0); // adds to history
    return 1;
}

int DataReducerMC::degreeTwoApplication(NetworKit::node node) {
    if (graph_.degree(node) != 2)
        return 0;
    return dominatingEdgeApplication(node); // adds to history
}

int DataReducerMC::dominatingEdgeApplication(NetworKit::node u) {
    assert(graph_.hasNode(u));
    auto degree = graph_.degree(u);
    assert(degree > 0);

    auto [highestAbsWeightNeighbor, highestAbsWeight] = graph_.getIthNeighborWithWeight(u, 0);
    NetworKit::edgeweight absSum = std::abs(highestAbsWeight);
    for (unsigned int i = 1; i < degree; i++) {
        auto [currentNeighbor, currentWeight] = graph_.getIthNeighborWithWeight(u, i);
        absSum += std::abs(currentWeight);
        if (std::abs(currentWeight) > std::abs(highestAbsWeight)) {
            highestAbsWeight = currentWeight;
            highestAbsWeightNeighbor = currentNeighbor;
        }
    }
    if (2 * std::abs(highestAbsWeight) >= absSum) {
        bool toSame = (highestAbsWeight <= 0);
        contractNodes(u, highestAbsWeightNeighbor, toSame); // adds to history
        numDominatingEdge_ += 1;
        return 1;
    }

    return 0;
}

int DataReducerMC::degreeThreeApplication(NetworKit::node node) {
    if (graph_.degree(node) != 3 || emphasis_ < 3) {
        return 0;
    }

    auto [neighborZero, weightZero] = graph_.getIthNeighborWithWeight(node, 0);
    auto [neighborOne, weightOne] = graph_.getIthNeighborWithWeight(node, 1);
    auto [neighborTwo, weightTwo] = graph_.getIthNeighborWithWeight(node, 2);

    std::array<NetworKit::node, 3> neighbors = {neighborZero, neighborOne, neighborTwo};

    // case 0: {0,1}, {0,2} and {1,2} are uncut, corresponds to coloring 0 0 0
    std::array<uint8_t, 4> cutCasesColors;
    std::array<NetworKit::edgeweight, 4> cutCasesValues;
    if (weightZero + weightOne + weightTwo > 0) {
        cutCasesColors[0] = 1;
        cutCasesValues[0] = weightZero + weightOne + weightTwo;
    } else {
        cutCasesColors[0] = 0;
        cutCasesValues[0] = 0;
    }
    // case 1: cut {0,1} and {0,2} uncut {1,2}, corresponds to coloring 0 1 1 or 1 0 0
    if (weightZero > weightOne + weightTwo) {
        cutCasesColors[3] = 1;
        cutCasesValues[1] = weightZero;
    } else {
        cutCasesColors[3] = 0;
        cutCasesValues[1] = weightOne + weightTwo;
    }
    // case 2: cut {0,1} and {1,2} uncut {0,2}, corresponds to coloring 0 1 0 or 1 0 1
    if (weightOne > weightZero + weightTwo) {
        cutCasesColors[2] = 0;
        cutCasesValues[2] = weightOne;
    } else {
        cutCasesColors[2] = 1;
        cutCasesValues[2] = weightZero + weightTwo;
    }
    // case 3: cut {0,2} and {1,2} uncut {0,1}, corresponds to coloring 0 0 1 or 1 1 0
    if (weightTwo > weightZero + weightOne) {
        cutCasesColors[1] = 0;
        cutCasesValues[3] = weightTwo;
    } else {
        cutCasesColors[1] = 1;
        cutCasesValues[3] = weightZero + weightOne;
    }

    NodeAggregation result = columnToNodeAggregation(node, neighbors, cutCasesColors);

    // update edge weights
    NetworKit::edgeweight oldWeights[3] = {0.0, 0.0, 0.0};
    if (graph_.hasEdge(neighborZero, neighborOne)) {
        oldWeights[0] = graph_.weight(neighborZero, neighborOne);
        graph_.removeEdge(neighborZero, neighborOne);
    }
    if (graph_.hasEdge(neighborZero, neighborTwo)) {
        oldWeights[1] = graph_.weight(neighborZero, neighborTwo);
        graph_.removeEdge(neighborZero, neighborTwo);
    }
    if (graph_.hasEdge(neighborOne, neighborTwo)) {
        oldWeights[2] = graph_.weight(neighborOne, neighborTwo);
        graph_.removeEdge(neighborOne, neighborTwo);
    }
    // perform updates
    auto newWeightZeroOne =
        oldWeights[0] + (cutCasesValues[1] + cutCasesValues[2] - cutCasesValues[3] - cutCasesValues[0]) / 2;
    if (newWeightZeroOne != 0)
        addEdge(neighborZero, neighborOne, newWeightZeroOne);
    auto newWeightZeroTwo =
        oldWeights[1] + (cutCasesValues[3] + cutCasesValues[1] - cutCasesValues[2] - cutCasesValues[0]) / 2;
    if (newWeightZeroTwo != 0)
        addEdge(neighborZero, neighborTwo, newWeightZeroTwo);
    auto newWeightOneTwo =
        oldWeights[2] + (cutCasesValues[3] + cutCasesValues[2] - cutCasesValues[1] - cutCasesValues[0]) / 2;
    if (newWeightOneTwo != 0)
        addEdge(neighborOne, neighborTwo, newWeightOneTwo);
    offset_ += cutCasesValues[0];

    removeVertex(node);
    for (auto curNeighbor : {neighborZero, neighborOne, neighborTwo}) {
        smartActivateAll(curNeighbor);
    }

    history_.push_back(result);
    degThreeReductions_ += 1;
    return 1;
}

// MightDo: Speed up, by early return
int DataReducerMC::threePathApplication(NetworKit::node u) {
    if (!allWeightsAreOne_ || graph_.degree(u) != 2)
        return 0; // no three path possible

    NetworKit::node neighbor0 = graph_.getIthNeighbor(u, 0);
    NetworKit::node neighbor1 = graph_.getIthNeighbor(u, 1);

    std::array<NetworKit::node, 4> path = {NetworKit::none, u, NetworKit::none, NetworKit::none};

    if (graph_.degree(neighbor0) == 2) {
        path[0] = neighbor1;
        path[2] = neighbor0;
        path[3] = u == graph_.getIthNeighbor(neighbor0, 1) ? graph_.getIthNeighbor(neighbor0, 0)
                                                           : graph_.getIthNeighbor(neighbor0, 1);
    }

    if (graph_.degree(neighbor1) == 2) {
        path[0] = neighbor0;
        path[2] = neighbor1;
        path[3] = u == graph_.getIthNeighbor(neighbor1, 1) ? graph_.getIthNeighbor(neighbor1, 0)
                                                           : graph_.getIthNeighbor(neighbor1, 1);
    }

    if ((path[0] != NetworKit::none) && (path[0] != path[3]) && (!graph_.hasEdge(path[0], path[3]))) {
        assert(!graph_.hasEdge(path[0], path[3]));
        removeVertex(path[1]);
        history_.emplace_back(path[1], path[0], true);
        removeVertex(path[2]);
        history_.emplace_back(path[2], path[3], true);

        addEdge(path[0], path[3], 1.0);
        smartActivateAll(path[0]);
        smartActivateAll(path[3]);
        // the new edge might make path[0] and path[3] external nodes of a clique
        // MightDo: The new clique would contain these two nodes, activating clique results in unnecessary checks
        // Still, it is enough, to activate the neighbors of one of the two vertices, as only shared ones are
        // relevant
        auto vertexWithSmallerDegree = graph_.degree(path[0]) < graph_.degree(path[3]) ? path[0] : path[3];
        for (auto neighbor : graph_.neighborRange(vertexWithSmallerDegree)) {
            activateRule(neighbor, RuleFlag::kWeightStableClique);
        }

        offset_ += 2;
        return 2;
    }

    return 0;
}

int DataReducerMC::cliqueApplication(NetworKit::node u) {
    assert(graph_.hasNode(u));
    if (!allWeightsAreOne_ || (graph_.degree(u) < 1))
        return 0;

    auto potentialCliqueSize = graph_.degree(u) + 1;
    auto potentialNearCliqueSize = graph_.degree(u) + 2;

    bool maybeInternalForClique = true;
    bool maybeInternalForNearClique = true;
    unsigned int numExternalForClique = 0;
    unsigned int numExternalForNearClique = 0;
    unsigned int numInternalForNearClique = 2; // the vertex and u and the vertex v for which {u,v} is missing
    std::vector<NetworKit::node> maybeInternalForCliqueVertices;
    maybeInternalForCliqueVertices.push_back(u);
    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        auto w = graph_.getIthNeighbor(u, i);
        if (graph_.degree(w) >= potentialCliqueSize) {
            if (graph_.degree(w) >= potentialNearCliqueSize)
                numExternalForNearClique += 1;
            else
                numInternalForNearClique += 1; // exactly one more neighbor
            numExternalForClique++;
            deactivateRule(w, RuleFlag::kWeightStableClique); // can not be internal node
        } else if (graph_.degree(w) == graph_.degree(u)) {
            deactivateRule(w, RuleFlag::kWeightStableClique); // it is enough to check current u
            maybeInternalForCliqueVertices.push_back(w);
            maybeInternalForNearClique = false;
        } else {
            // If u has a neighbor with lower degree its definitely not an internal vertex of a clique
            maybeInternalForClique = false;
            maybeInternalForNearClique = false;
        }
    }

    // start with clique to near clique first, is faster
    if (maybeInternalForClique                                      // reduction might be possible
        && (numExternalForClique > divTwoCeil(potentialCliqueSize)) // to many external for clique reduction
        && (maybeInternalForCliqueVertices.size() > 2)              // for two internal, twin-rule will catch this later
        && neighborhoodIsClique(u, true)) {
        graph_.removeEdge(maybeInternalForCliqueVertices[0], maybeInternalForCliqueVertices[1]);
        // offset does not change
        return 1;
    }

    // real clique check
    if (maybeInternalForClique                                       // might be internal for clique
        && (numExternalForClique <= divTwoCeil(potentialCliqueSize)) // few enough external
        && neighborhoodIsClique(u, true)) {
        numCliques_ += 1;
        return cliqueReductionCore(u, potentialCliqueSize, numExternalForClique);
    }

    // check for near-clique which allows for clique reduction via implicit adding of edge to form full clique
    if (maybeInternalForNearClique
        && (numExternalForNearClique <= divTwoCeil(potentialNearCliqueSize))      // reduction might be possible
        && ((potentialNearCliqueSize % 2 == 1) || (numInternalForNearClique > 2)) // additional restriction
        && neighborhoodIsClique(u, true)) {
        // here, all neighbors of u are guaranteed to have higher degree
        auto neighbor = graph_.getIthNeighbor(u, 0); // neighbor with smallest degree would be better
        for (auto v : graph_.neighborRange(u)) {
            neighborWeights_[v] = 1.0;
        }
        auto partnerVertex = NetworKit::none;
        for (auto v : graph_.neighborRange(neighbor)) {
            if (u != v && (graph_.degree(v) == graph_.degree(u))) { // TODO: check if clique is active?
                // u and v might be internal for a near clique
                bool sameNeighborhood = true;
                for (auto w : graph_.neighborRange(v)) {
                    sameNeighborhood &= (neighborWeights_[w] == 1.0);
                }
                if (sameNeighborhood) {
                    partnerVertex = v;
                    break;
                }
            }
        }
        // cleanup
        for (auto v : graph_.neighborRange(u)) {
            neighborWeights_[v] = 0.0;
        }
        if (partnerVertex != NetworKit::none) {
            numNearCliques_ += 1;
            graph_.addEdge(u, partnerVertex, 1.0); // all weights are one
            return nearCliqueReductionCore(u, partnerVertex, potentialNearCliqueSize, numExternalForNearClique);
        }
    }

    return 0;
}

int DataReducerMC::cliqueReductionCore(NetworKit::node u, unsigned int cliqueSize, unsigned int numExternalNodes) {
    const auto numInternalNodes = cliqueSize - numExternalNodes;
    std::vector<NetworKit::node> internalNodes = {u};
    internalNodes.reserve(numInternalNodes);
    auto externalNodes = std::vector<NetworKit::node>();
    externalNodes.reserve(numExternalNodes);

    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        auto w = graph_.getIthNeighbor(u, i);
        if (graph_.degree(w) > graph_.degree(u)) {
            externalNodes.push_back(w);
        } else if (graph_.degree(w) == graph_.degree(u)) {
            internalNodes.push_back(w);
        }
    }
    assert(internalNodes.size() == numInternalNodes);
    assert(externalNodes.size() == numExternalNodes);

    if (numExternalNodes == 0) {
        for (unsigned int i = 0; i < internalNodes.size() - 1; i++) {
            history_.emplace_back(internalNodes[i], internalNodes[i + 1], true);
            removeVertex(internalNodes[i]);
        }
        degreeZeroApplication(internalNodes.back());
        assert(std::all_of(internalNodes.begin(), internalNodes.end(), [&](auto u) { return !graph_.hasNode(u); }));
        auto sizeSplit = cutSizePerfectSplit(cliqueSize);
        offset_ += sizeSplit;
        return static_cast<int>(numInternalNodes);
    }

    for (unsigned int i = numInternalNodes; i > numExternalNodes; i--) {
        history_.emplace_back(internalNodes[i - 1], internalNodes[i - 2], true);
        removeVertex(internalNodes[i - 1]);
    }
    for (unsigned int i = 0;
         i < std::min(static_cast<uint64_t>(numInternalNodes), static_cast<uint64_t>(numExternalNodes)); i++) {
        history_.emplace_back(internalNodes[i], externalNodes[i], true);
        removeVertex(internalNodes[i]);
    }

    for (unsigned int i = 0; i < externalNodes.size(); i++) {
        auto ext1 = externalNodes[i];
        for (unsigned int j = i + 1; j < externalNodes.size(); j++) {
            auto ext2 = externalNodes[j];
            assert(graph_.hasEdge(ext1, ext2));
            graph_.removeEdge(ext1, ext2);
        }
        for (auto neighbor : graph_.neighborRange(ext1)) {
            if (graph_.degree(neighbor) == 2)
                activateRule(neighbor, RuleFlag::kWeightStableClique);
        }
    }

    for (auto externalNode : externalNodes) {
        smartActivateAll(externalNode); // separate loop here, as neighborhood needs to be final
    }

    offset_ += cutSizePerfectSplit(cliqueSize);
    return static_cast<int>(numInternalNodes);
}

int DataReducerMC::nearCliqueReductionCore(NetworKit::node u, NetworKit::node partner, unsigned int cliqueSize,
                                           unsigned int numExternalNodes) {
    const auto numInternalNodes = cliqueSize - numExternalNodes;
    assert(numInternalNodes + 1 >= numExternalNodes);
    assert(cliqueSize >= 3);

    auto internalNodes = std::vector<NetworKit::node>();
    internalNodes.reserve(numInternalNodes);
    auto externalNodes = std::vector<NetworKit::node>();
    externalNodes.reserve(numExternalNodes);

    internalNodes.push_back(u);
    internalNodes.push_back(partner);
    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        auto w = graph_.getIthNeighbor(u, i);
        if (graph_.degree(w) > graph_.degree(u)) {
            externalNodes.push_back(w);
        } else if ((graph_.degree(w) == graph_.degree(u)) && (w != partner)) {
            internalNodes.push_back(w);
        }
    }
    assert(internalNodes.size() == numInternalNodes);
    assert(externalNodes.size() == numExternalNodes);
    assert(numInternalNodes + numExternalNodes == cliqueSize);

    auto cliqueVerticesLeft = externalNodes;
    assert(cliqueVerticesLeft.size() == externalNodes.size());
    for (auto intVertex : internalNodes) {
        cliqueVerticesLeft.push_back(intVertex);
    }
    assert(cliqueVerticesLeft.size() == cliqueSize);
    assert(cliqueVerticesLeft[numExternalNodes] == u);
    assert(cliqueVerticesLeft[numExternalNodes + 1] == partner);
    for (int i = numInternalNodes - 1; i >= 2; i--) {
        cliqueVerticesLeft.pop_back();
        // to smaller of all before
        history_.emplace_back(internalNodes[i], cliqueVerticesLeft);
        removeVertex(internalNodes[i]);
    }
    // partner to same
    history_.emplace_back(partner, u, false);
    removeVertex(partner);
    // u goes to smaller of all external
    history_.emplace_back(u, externalNodes);
    removeVertex(u);

    for (unsigned int i = 0; i < externalNodes.size(); i++) {
        auto ext1 = externalNodes[i];
        for (unsigned int j = i + 1; j < externalNodes.size(); j++) {
            auto ext2 = externalNodes[j];
            assert(graph_.hasEdge(ext1, ext2));
            graph_.removeEdge(ext1, ext2);
        }
        for (auto neighbor : graph_.neighborRange(ext1)) {
            if (graph_.degree(neighbor) == 2)
                activateRule(neighbor, RuleFlag::kWeightStableThreePath);
        }
    }

    for (auto externalNode : externalNodes) {
        smartActivateAll(externalNode); // separate loop here, as neighborhood needs to be final
    }

    offset_ += cutSizePerfectSplit(cliqueSize);
    return static_cast<int>(numInternalNodes);
}

int DataReducerMC::gomoryHuApplication() {
    if (graph_.numberOfNodes() < 3)
        return 0;

    // copy of graph with absolute edge weights
    auto g = graph_;
    g.indexEdges(true);

    for (auto edge : graph_.edgeWeightRange()) {
        if (edge.weight < 0)
            g.setWeight(edge.u, edge.v, std::abs(edge.weight));
    }

    int numContractions = 0;
    auto contractions = std::vector<std::vector<NetworKit::node>>(graph_.upperNodeIdBound());
    auto contractedToNode = std::vector<NetworKit::node>(graph_.upperNodeIdBound());
    for (auto node : g.nodeRange()) {
        contractedToNode[node] = node;
        contractions[node] = {node};
    }

    auto ght = GomoryHuTree(g);
    ght.run();

    for (auto edge : g.edgeWeightRange()) {
        auto ghtEdge = ght.minCutEdge(edge.u, edge.v);
        // check if edge is dominating
        if (2 * edge.weight >= ghtEdge.weight) {
            auto nodeToDelete = contractedToNode[edge.u];
            auto nodeToKeep = contractedToNode[edge.v];
            if (nodeToDelete == nodeToKeep)
                continue;
            bool toSame = (graph_.weight(nodeToDelete, nodeToKeep) <= 0);
            contractNodes(nodeToDelete, nodeToKeep, toSame); // adds to history
            numContractions++;
            for (auto node : contractions[nodeToDelete]) {
                contractedToNode[node] = contractedToNode[edge.v];
                contractions[nodeToKeep].push_back(node);
            }
            contractions[nodeToDelete] = {};
            // contract edge in gomory hu tree
            ght.contractNodesInGHT(ghtEdge.u, ghtEdge.v);
        }
    }
    numGomoryHuReductions_ += numContractions;
    return numContractions;
}

int DataReducerMC::smallSeparatorApplication() {
    auto timer = Aux::Timer();
    timer.start();
    int removedNodes = 0;
    auto separator = FastVertexSeparator(graph_);
    separator.run(20);
    if (!separator.separatedSets().empty()) {
        for (unsigned int i = 0; i < separator.numberOfSeparators(); i++) {
            const std::vector<NetworKit::node> &separatedSet = separator.separatedSets()[i];
            const std::vector<NetworKit::node> &sep = separator.separators()[i];
            if (!std::all_of(separatedSet.begin(), separatedSet.end(), [&](auto v) { return graph_.hasNode(v); })
                || !std::all_of(sep.begin(), sep.end(), [&](auto v) { return graph_.hasNode(v); }))
                continue;
            switch (sep.size()) {
                case 0:
                    break;
                case 1:
                    removedNodes += oneSepApplication(separatedSet, sep);
                    break;
                case 2:
                    removedNodes += twoSepApplication(separatedSet, sep);
                    break;
                case 3:
                    if (emphasis_ > 2)
                        removedNodes += threeSepApplication(separatedSet, sep);
                    break;
                default:
                    assert(false);
            }
        }
    }
    numSmallSeparatorReductions_ += removedNodes;
    return removedNodes;
}

int DataReducerMC::oneSepApplication(const std::vector<NetworKit::node> &separatedSet,
                                     const std::vector<NetworKit::node> &sep) {
    std::vector<NetworKit::node> vs = sep;
    vs.insert(vs.end(), separatedSet.begin(), separatedSet.end());
    auto cp = inducedSubgraphCompact(graph_, vs);
    // calculate two solutions
    auto es = EnumerationSolver(cp.compactGraph);
    es.run({cp.orig2compact.at(sep[0])}, {false});
    offset_ += es.bestSolutionValue();

    // update history
    for (auto j : separatedSet) {
        auto cpVertexId = cp.orig2compact.at(j);
        history_.emplace_back(j, sep[0], es.partitionOf(cpVertexId) == 1);
        removeVertex(j);
    }
    for (auto v : sep)
        smartActivateAll(v);
    return static_cast<int>(separatedSet.size());
}

int DataReducerMC::twoSepApplication(const std::vector<NetworKit::node> &separatedSet,
                                     const std::vector<NetworKit::node> &sep) {
    std::vector<NetworKit::node> vs = sep;
    vs.insert(vs.end(), separatedSet.begin(), separatedSet.end());
    auto cp = inducedSubgraphCompact(graph_, vs);
    // calculate two solutions
    auto esSame = EnumerationSolver(cp.compactGraph);
    esSame.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1])}, {false, false});
    auto esDiff = EnumerationSolver(cp.compactGraph);
    esDiff.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1])}, {false, true});
    // update edge weights
    if (graph_.hasEdge(sep[0], sep[1])) {
        graph_.removeEdge(sep[0], sep[1]);
    }
    offset_ += esSame.bestSolutionValue();
    auto newEdgeWeight = esDiff.bestSolutionValue() - esSame.bestSolutionValue();
    if (newEdgeWeight != 0) {
        addEdge(sep[0], sep[1], newEdgeWeight);
    }

    // update history
    for (auto j : separatedSet) {
        auto cpVertexId = cp.orig2compact.at(j);
        std::array<uint8_t, 2> col = {esSame.partitionOf(cpVertexId), esDiff.partitionOf(cpVertexId)};
        history_.push_back(columnToNodeAggregation(j, {sep[0], sep[1]}, col));
        removeVertex(j);
    }
    for (auto v : sep)
        smartActivateAll(v);

    return static_cast<int>(separatedSet.size());
}

int DataReducerMC::threeSepApplication(const std::vector<NetworKit::node> &separatedSet,
                                       const std::vector<NetworKit::node> &sep) {
    // get induced subgraph
    std::vector<NetworKit::node> vs = sep;
    vs.insert(vs.end(), separatedSet.begin(), separatedSet.end());
    auto cp = inducedSubgraphCompact(graph_, vs);

    auto es000 = EnumerationSolver(cp.compactGraph);
    es000.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1]), cp.orig2compact.at(sep[2])},
              {false, false, false});

    auto es001 = EnumerationSolver(cp.compactGraph);
    es001.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1]), cp.orig2compact.at(sep[2])},
              {false, false, true});

    auto es010 = EnumerationSolver(cp.compactGraph);
    es010.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1]), cp.orig2compact.at(sep[2])},
              {false, true, false});

    auto es011 = EnumerationSolver(cp.compactGraph);
    es011.run({cp.orig2compact.at(sep[0]), cp.orig2compact.at(sep[1]), cp.orig2compact.at(sep[2])},
              {false, true, true});

    offset_ += es000.bestSolutionValue();

    // update edge weights
    if (graph_.hasEdge(sep[0], sep[1])) {
        graph_.removeEdge(sep[0], sep[1]);
    }

    if (graph_.hasEdge(sep[0], sep[2])) {
        graph_.removeEdge(sep[0], sep[2]);
    }

    if (graph_.hasEdge(sep[1], sep[2])) {
        graph_.removeEdge(sep[1], sep[2]);
    }
    std::array<NetworKit::edgeweight, 4> cutCasesValues = {es000.bestSolutionValue(), es011.bestSolutionValue(),
                                                           es010.bestSolutionValue(), es001.bestSolutionValue()};
    auto newWeightZeroOne = (cutCasesValues[1] + cutCasesValues[2] - cutCasesValues[3] - cutCasesValues[0]) / 2;
    if (newWeightZeroOne != 0)
        addEdge(sep[0], sep[1], newWeightZeroOne);
    auto newWeightZeroTwo = (cutCasesValues[3] + cutCasesValues[1] - cutCasesValues[2] - cutCasesValues[0]) / 2;
    if (newWeightZeroTwo != 0)
        addEdge(sep[0], sep[2], newWeightZeroTwo);
    auto newWeightOneTwo = (cutCasesValues[3] + cutCasesValues[2] - cutCasesValues[1] - cutCasesValues[0]) / 2;
    if (newWeightOneTwo != 0)
        addEdge(sep[1], sep[2], newWeightOneTwo);

    // update history
    for (auto j : separatedSet) {
        auto cpVertexId = cp.orig2compact.at(j);
        std::array<uint8_t, 4> col = {es000.partitionOf(cpVertexId), es001.partitionOf(cpVertexId),
                                      es010.partitionOf(cpVertexId), es011.partitionOf(cpVertexId)};
        history_.push_back(columnToNodeAggregation(j, {sep[0], sep[1], sep[2]}, col));
        removeVertex(j);
    }
    for (auto v : sep)
        smartActivateAll(v);
    return static_cast<int>(separatedSet.size());
}

void DataReducerMC::contractNodes(NetworKit::node nodeToDelete, NetworKit::node nodeToKeep, bool toSamePartition) {
    assert(nodeToKeep != nodeToDelete);
    assert(graph_.hasNode(nodeToDelete));
    assert(graph_.hasNode(nodeToKeep));
    assert(std::all_of(neighborWeights_.begin(), neighborWeights_.end(), [](auto w) { return w == 0.0; }));

    if (!toSamePartition)
        offset_ += graph_.weightedDegree(nodeToDelete);

    for (auto [currentNeighbor, currentWeight] : graph_.weightNeighborRange(nodeToKeep)) {
        neighborWeights_[currentNeighbor] = currentWeight;
    }

    for (auto [currentNeighbor, currentWeight] : graph_.weightNeighborRange(nodeToDelete)) {
        assert(graph_.hasNode(currentNeighbor));
        assert(currentNeighbor != nodeToDelete); // no self loop

        if (currentNeighbor != nodeToKeep) {
            NetworKit::edgeweight oldEdgeWeight = neighborWeights_[currentNeighbor];
            neighborWeights_[currentNeighbor] = 0.0; // reset here as edge might get deleted

            if (oldEdgeWeight != 0.0) {
                graph_.removeEdge(currentNeighbor, nodeToKeep);
            }

            if (toSamePartition && ((oldEdgeWeight + currentWeight) != 0)) {
                allWeightsAreOne_ &= (oldEdgeWeight + currentWeight == 1.);
                addEdge(nodeToKeep, currentNeighbor, (oldEdgeWeight + currentWeight));
            } else if (!toSamePartition && ((oldEdgeWeight - currentWeight) != 0)) {
                allWeightsAreOne_ &= (oldEdgeWeight - currentWeight == 1.);
                addEdge(nodeToKeep, currentNeighbor, (oldEdgeWeight - currentWeight));
            }

            activateRule(currentNeighbor, RuleFlag::kAll); // MightDo: Can we use smart activate here?
        }
    }

    // final cleanup
    for (auto currentNeighbor : graph_.neighborRange(nodeToKeep)) {
        neighborWeights_[currentNeighbor] = 0.0;
    }
    assert(std::all_of(neighborWeights_.begin(), neighborWeights_.end(), [](auto w) { return w == 0.0; }));

    deactivateAllRules(nodeToDelete);
    history_.emplace_back(nodeToDelete, nodeToKeep, !toSamePartition);
    removeVertex(nodeToDelete);
    smartActivateAll(nodeToKeep);
}

/*
 * Given a full solution for the result of the data reduction, recover a solution for the original graph
 */
std::vector<Partition> DataReducerMC::recoverSolution(const std::vector<Partition> &transformedSol) const {
    return sms::recoverSolution(history_, transformedSol);
}

int DataReducerMC::triangleApplication(NetworKit::node u) {
    assert(std::all_of(neighborMarks_.begin(), neighborMarks_.end(), [](auto x) { return x == NetworKit::none; }));
    if (allWeightsAreOne_ || !useTriangles_)
        return 0;

    auto neighborRangeU = graph_.neighborRange(u);
    auto neighborsU = std::vector<NetworKit::node>(neighborRangeU.begin(), neighborRangeU.end());
    for (size_t i = 0; i < graph_.degree(u); i++)
        neighborMarks_[graph_.getIthNeighbor(u, i)] = i;

    int cnt = 0;
    for (auto [v, uvWeight] : graph_.weightNeighborRange(u)) {
        neighborMarks_[v] = NetworKit::none; // reset mark
        if (isRuleActive(v, RuleFlag::kTriangle))
            continue;

        for (auto [w, vwWeight] : graph_.weightNeighborRange(v)) {
            if ((w == u) || isRuleActive(w, RuleFlag::kTriangle)
                || (weightOneNeighborhood_[u] && weightOneNeighborhood_[v] && weightOneNeighborhood_[w])
                || (neighborMarks_[w] == NetworKit::none))
                continue;
            assert(graph_.getIthNeighbor(u, neighborMarks_[w]) == w);
            auto uwWeight = graph_.getIthNeighborWeight(u, neighborMarks_[w]);
            auto weightedTriangle = WeightedTriangle{
                .u = u, .v = v, .w = w, .uvWeight = uvWeight, .vwWeight = vwWeight, .wuWeight = uwWeight};
            auto res = triangleApplicationCore(weightedTriangle);
            cnt += 1;
            if (res > 0 || cnt > 20) // ToDo: No hardcoded limit for triangles
            {
                for (auto x : neighborsU)
                    neighborMarks_[x] = NetworKit::none; // reset early, because we return here
                numTriangles_ += 1;
                return res;
            }
        }
    }

    return 0;
}

int DataReducerMC::triangleApplicationCore(WeightedTriangle &wt) {
    assert(graph_.hasEdge(wt.u, wt.v) && graph_.hasEdge(wt.v, wt.w) && graph_.hasEdge(wt.w, wt.u));

    int numberOfNegativeEdges = (wt.uvWeight < 0) + (wt.vwWeight < 0) + (wt.wuWeight < 0);

    switch (numberOfNegativeEdges) {
        case 3: {
            NetworKit::edgeweight absOutU = absOutWeight(graph_, wt.u);
            NetworKit::edgeweight absOutV = absOutWeight(graph_, wt.v);
            NetworKit::edgeweight absOutW = absOutWeight(graph_, wt.w);

            if (canContractThreeNegative(wt.uvWeight, wt.wuWeight, wt.vwWeight, absOutU, absOutV, absOutW)) {
                contractNodes(wt.u, wt.v, true);
                return 1;
            } else if (canContractThreeNegative(wt.wuWeight, wt.vwWeight, wt.uvWeight, absOutW, absOutU, absOutV)) {
                contractNodes(wt.u, wt.w, true);
                return 1;
            } else if (canContractThreeNegative(wt.vwWeight, wt.uvWeight, wt.wuWeight, absOutV, absOutW, absOutU)) {
                contractNodes(wt.v, wt.w, true);
                return 1;
            }
            return 0;
        }
        case 1: {
            // make sure edge {u, v} is the negative edge
            if (wt.wuWeight < 0) {
                std::swap(wt.uvWeight, wt.wuWeight);
                std::swap(wt.v, wt.w);
            } else if (wt.vwWeight < 0) {
                std::swap(wt.uvWeight, wt.vwWeight);
                std::swap(wt.u, wt.w);
            }
            assert(wt.uvWeight < 0);

            NetworKit::edgeweight absOutU = absOutWeight(graph_, wt.u);
            NetworKit::edgeweight absOutV = absOutWeight(graph_, wt.v);
            NetworKit::edgeweight absOutW = absOutWeight(graph_, wt.w);

            if (useTriangleThree_
                && canContractThreeNegative(wt.uvWeight, -wt.wuWeight, -wt.vwWeight, absOutU, absOutV, absOutW)) {
                contractNodes(wt.u, wt.v, true);
                // numTriangleThreeContractions_ += 1;
                return 1;
            } else if (canContractThreeNegative(-wt.wuWeight, -wt.vwWeight, wt.uvWeight, absOutW, absOutU, absOutV)) {
                contractNodes(wt.u, wt.w, false);
                // numTriangleTwoContractions_ += 1;
                return 1;
            } else if (canContractThreeNegative(-wt.vwWeight, wt.uvWeight, -wt.wuWeight, absOutV, absOutW, absOutU)) {
                contractNodes(wt.v, wt.w, false);
                // numTriangleTwoContractions_ += 1;
                return 1;
            }
            return 0;
        }
        default:
            return 0;
    }
}

NodeAggregation columnToNodeAggregation(NetworKit::node u, const std::array<NetworKit::node, 2> &separator,
                                        const std::array<uint8_t, 2> &column) {
    assert(std::all_of(column.begin(), column.end(), [](auto x) { return x == 0 || x == 1; }));

    int numberOfOnes = column[0] + column[1];
    NetworKit::node aggregateTo = separator[numberOfOnes % 2];
    bool same = (column[0] == 0);
    return NodeAggregation(u, aggregateTo, !same);
}

NodeAggregation columnToNodeAggregation(const NetworKit::node u, const std::array<NetworKit::node, 3> &separator,
                                        const std::array<uint8_t, 4> &column) {
    assert(std::all_of(column.begin(), column.end(), [](auto x) { return x == 0 || x == 1; }));

    auto numberOfOnes = std::count(column.begin(), column.end(), 1);

    switch (numberOfOnes) {
        case 0:
            return NodeAggregation(u, separator[0], false);
        case 1:
            if (column[0] == 1) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], false, true, true);
            } else if (column[1] == 1) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], false, true, false);
            } else if (column[2] == 1) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], false, false, true);
            } else {
                assert(column[3] == 1);
                return NodeAggregation(u, separator[0], separator[1], separator[2], false, false, false);
            }
        case 2:
            if (column[0] == 0) {
                if (column[1] == 0) { // 00
                    return NodeAggregation(u, separator[1], false);
                } else { // 01
                    assert(column[2] == 0);
                    return NodeAggregation(u, separator[2], false);
                }
            } else {
                assert(column[0] == 1);
                if (column[1] == 1) { // 11
                    return NodeAggregation(u, separator[1], true);
                } else { // 10
                    assert(column[2] == 1);
                    return NodeAggregation(u, separator[2], true);
                }
            }
        case 3:
            if (column[0] == 0) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], true, true, true);
            } else if (column[1] == 0) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], true, true, false);
            } else if (column[2] == 0) {
                return NodeAggregation(u, separator[0], separator[1], separator[2], true, false, true);
            } else {
                assert(column[3] == 0);
                return NodeAggregation(u, separator[0], separator[1], separator[2], true, false, false);
            }
        case 4:
            return NodeAggregation(u, separator[0], true);
        default:
            throw std::runtime_error("Invalid number of ones in column");
    }
}

std::vector<Partition> recoverSolution(const std::vector<NodeAggregation> &aggr,
                                       const std::vector<Partition> &transformedSol) {
    std::vector<Partition> res(transformedSol);
    for (unsigned int i = aggr.size(); i > 0; i--) { // in reverse order of NodeAggregations
        const auto &curAggr = aggr[i - 1];
        auto nodeToPartition = curAggr.toAggregate;
        assert(nodeToPartition < transformedSol.size());
        assert(transformedSol[nodeToPartition] == Partition::kUNASSIGNED);
        assert(res[nodeToPartition] == Partition::kUNASSIGNED);
        Partition recoveredPartition;
        std::vector<bool> relevantPartitions(curAggr.others.size(), false);
        for (int j = 0; j < std::ssize(curAggr.others); j++) {
            assert(res[curAggr.others[j]] != Partition::kUNASSIGNED);
            relevantPartitions[j] = static_cast<bool>(res[curAggr.others[j]]);
        }
        recoveredPartition = static_cast<Partition>(curAggr.reverse(relevantPartitions));
        res[nodeToPartition] = recoveredPartition;
    }
    return res;
}

bool DataReducerMC::neighborhoodIsClique(NetworKit::node v, bool degreesChecked) {
    assert(graph_.hasNode(v));
    assert(graph_.hasNode(v));
    assert(std::all_of(neighborWeights_.begin(), neighborWeights_.end(), [](auto i) { return i == 0.0; }));

    if (not degreesChecked) {
        for (NetworKit::node w : graph_.neighborRange(v)) {
            if (graph_.degree(w) < graph_.degree(v))
                return false;
        }
    }

    bool isGood = true;
    unsigned int cliqueSize = graph_.degree(v) + 1;

    neighborWeights_[v] = 1.0;
    for (NetworKit::node w : graph_.neighborRange(v))
        neighborWeights_[w] = 1.0;

    for (NetworKit::node w : graph_.neighborRange(v)) {
        NetworKit::edgeweight numSharedNeighbors = 0;
        // Count the number of clique neighbors of w
        for (NetworKit::node x : graph_.neighborRange(w)) {
            numSharedNeighbors += neighborWeights_[x];
        }
        // Check if w neighbors every node in the supposed clique
        if (numSharedNeighbors != cliqueSize - 1) {
            isGood = false;
            break;
        }
    }

    // Restore the markings
    neighborWeights_[v] = 0.0;
    for (NetworKit::node w : graph_.neighborRange(v))
        neighborWeights_[w] = 0;
    assert(std::all_of(neighborWeights_.begin(), neighborWeights_.end(), [](auto i) { return i == 0.0; }));
    return isGood;
}

void DataReducerMC::smartActivateAll(NetworKit::node u) {
    RuleFlag rules = getDegreeFlag(u);
    rules = rules | getNonDegreeFlags(u);
    vertexStatus_[u] = vertexStatus_[u] | rules;
    vertexHeap_.update(u);
}

RuleFlag DataReducerMC::getDegreeFlag(NetworKit::node u) const {
    auto rules = RuleFlag::kNONE;
    switch (graph_.degree(u)) {
        case 0:
            rules = rules | RuleFlag::kDegreeZero;
            break;
        case 1:
            rules = rules | RuleFlag::kDegreeOne;
            break;
        case 2:
            rules = rules | RuleFlag::kDegreeTwo;
            break;
        case 3:
            if (useDegreeThree_) {
                rules = rules | RuleFlag::kDegreeThree;
            }
            break;
        default:; // do nothing
    }
    return rules;
}

RuleFlag DataReducerMC::getNonDegreeFlags(NetworKit::node u) const {
    auto rules = RuleFlag::kNONE;
    if (weightOneNeighborhood_[u]) {
        rules = rules | RuleFlag::kWeightStableThreePath;
        rules = rules | RuleFlag::kWeightStableClique;
    } else {
        rules = rules | RuleFlag::kDominatingEdge;
        rules = rules | RuleFlag::kTriangle;
    }
    rules = rules | RuleFlag::kSimilarNeighborhood;
    return rules;
}

int DataReducerMC::similarNeighborhoodApplication(NetworKit::node u) {
    assert(std::all_of(neighborWeights2_.begin(), neighborWeights2_.end(), [](auto w) { return w == 0.; }));
    int res = 0;

    auto neighborRange = graph_.neighborRange(u);
    auto neighborsU = std::vector<NetworKit::node>(neighborRange.begin(), neighborRange.end());
    for (auto [v, weight] : graph_.weightNeighborRange(u)) {
        neighborWeights2_[v] = weight;
    }

    for (auto [v, weight] : graph_.weightNeighborRange(u)) {
        // the check will be performed the other way around later
        if (!isRuleActive(v, RuleFlag::kSimilarNeighborhood)) {
            auto alpha = neighborhoodAlphaMarksSet(graph_, u, v, true, neighborWeights2_);
            if (alpha > 0.0 && weight < 0.0) {
                posAlphaWithEdge_ += 1;
                contractNodes(v, u, true);
                res = 1;
                break;
            } else if (alpha < 0.0 && weight > 0.0) {
                negAlphaWithEdge_ += 1;
                contractNodes(v, u, false);
                res = 1;
                break;
            } else if (useTwins_ && alpha == 1.0 && weightOneNeighborhood_[u] && graph_.degree(u) % 2 == 0) {
                twins_ += 1;
                contractNodes(v, u, true);
                res = 1;
                break;
            }
        }
    }

    for (auto v : neighborsU)
        neighborWeights2_[v] = 0.0;

    return res;
}

int DataReducerMC::removeSeparatedCliques() {
    // hash vertices
    // key = sorted {neighborhood \setunion {themselves}}
    // value = list of vertices with this neighborhood
    auto buckets = std::map<std::vector<NetworKit::node>, std::vector<NetworKit::node>>();

    for (auto u : graph_.nodeRange()) {
        auto neighborRange = graph_.neighborRange(u);
        std::vector<NetworKit::node> neighbors(neighborRange.begin(), neighborRange.end());
        neighbors.push_back(u);
        std::sort(neighbors.begin(), neighbors.end());

        if (buckets.find(neighbors) == buckets.end()) {
            buckets[neighbors] = std::vector<NetworKit::node>{u};
        } else {
            buckets[neighbors].push_back(u);
        }
    }

    std::vector<bool> markings_(graph_.upperNodeIdBound(), 0);
    int numVerticesRemoved = 0;
    std::vector<bool> touched(graph_.upperNodeIdBound(), false);
    // loop over each bucket
    // for each one it is guaranteed that the nodes are in the same clique and have the same neighbors
    for (const auto &[cliqueAndNeighborhood, clique] : buckets) {
        if (!std::ranges::all_of(cliqueAndNeighborhood, [&touched](auto u) { return !touched[u]; }))
            continue;

        // check total size of neighborhood, which for any node u is deg(u) +1 - bucketSize
        // if neighborhood size is less or equal to bucket size, then the clique is removable
        if ((clique.size() > 1) && (2 * clique.size() + 1 >= cliqueAndNeighborhood.size())) {
            for (auto u : cliqueAndNeighborhood) {
                touched[u] = true;
            }

            // collect neighborhood of clique
            // mark nodes which are part of the clique
            for (auto u : clique) {
                markings_[u] = true;
            }
            // unmarked nodes form neighborhood
            std::vector<NetworKit::node> neighborhood;
            for (auto u : cliqueAndNeighborhood) {
                if (!markings_[u]) {
                    neighborhood.push_back(u);
                }
            }
            // reset marks
            for (auto u : clique) {
                markings_[u] = false;
            }

            if (useFullSeparatedCliques_) {
                for (unsigned int i = 0; i < neighborhood.size(); i++) {
                    for (unsigned int j = i + 1; j < neighborhood.size(); j++) {
                        if (graph_.hasEdge(neighborhood[i], neighborhood[j])) {
                            graph_.removeEdge(neighborhood[i], neighborhood[j]);
                        } else {
                            graph_.addEdge(neighborhood[i], neighborhood[j], -1.);
                        }
                    }
                }

                for (unsigned int i = clique.size(); i > neighborhood.size(); i--) {
                    history_.emplace_back(clique[i - 1], clique[i - 2], true);
                    removeVertex(clique[i - 1]);
                }
                for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(clique.size()),
                                                      static_cast<unsigned int>(neighborhood.size()));
                     i++) {
                    history_.emplace_back(clique[i], neighborhood[i], true);
                    removeVertex(clique[i]);
                }
                numVerticesRemoved += std::ssize(clique);
                offset_ += cutSizePerfectSplit(cliqueAndNeighborhood.size());
            } else {
                std::vector<NetworKit::node> remainingCliqueNodes = clique;
                const auto curNeighborhoodSize = cliqueAndNeighborhood.size() - remainingCliqueNodes.size();

                while (remainingCliqueNodes.size() > curNeighborhoodSize && (!remainingCliqueNodes.empty())) {
                    // remove two nodes from clique
                    auto x1 = remainingCliqueNodes.back();
                    remainingCliqueNodes.pop_back();
                    auto x2 = remainingCliqueNodes.back();
                    remainingCliqueNodes.pop_back();
                    contractNodes(x1, x2, false);
                    assert(graph_.hasNode(x2));
                    assert(graph_.degree(x2) == 0);
                    removeVertex(x2);
                    history_.emplace_back(x2);
                    numVerticesRemoved += 2;
                }
            }
        }
    }
    numSeparatedCliqueVerticesRemoved_ += numVerticesRemoved;
    return numVerticesRemoved;
}

} // namespace sms