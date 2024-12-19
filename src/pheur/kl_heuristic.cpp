#include "sms/pheur/kl_heuristic.hpp"

#include <chrono>

#include "tlx/container/d_ary_addressable_int_heap.hpp"

#include "sms/auxiliary/math.hpp"
#include "sms/graph/graphs.hpp"

namespace sms {

const std::vector<bool> &KLHeuristic::getPrimalSolution() const {
    return solution_;
}

void KLHeuristic::phase1Optimization(const std::vector<bool> &initSol,
                                     const std::vector<NetworKit::node> &fixedVertices) {
    assert(fixedVertices.size() <= originalGraph_->numberOfNodes());
    t_.start();

    fixVertices(fixedVertices);
    auto workingPartition = BiPartition(initSol); // copy to not change the input
    assert(workingPartition.valid());

    phase1OptimizationSubRoutine(workingPartition);

    for (size_t i = 0; i < originalGraph_->upperNodeIdBound(); i++) {
        solution_[i] = workingPartition.inA[i];
    }

    t_.stop();
}

void KLHeuristic::phase1OptimizationSubRoutine(BiPartition &partition) {
    bool improvement = true;
    int numOuterIterations = 0;

    while (improvement && belowMaxRuntime() && belowDualBound(partition.inA)) {
        std::vector<Swap> swapSuggestions;
        swapSuggestions.reserve(originalGraph_->numberOfNodes());
        auto numInnerIterations = 0;
        auto improvementValues = improvementValuesFromScratch(partition);

        const auto heapCompare = [&improvementValues](auto a, auto b) {
            return improvementValues[a] > improvementValues[b];
        };
        tlx::DAryAddressableIntHeap<NetworKit::node, 2, decltype(heapCompare)> heap(heapCompare);
        for (auto u : originalGraph_->nodeRange()) {
            if (!vertexIsFixed_[u])
                heap.push(u);
        }
        while ((!heap.empty()) && belowMaxRuntime()) {
            assert(swapSuggestions.size() == static_cast<size_t>(numInnerIterations));
            assert(heap.sanity_check());
            auto nextBestVertex = heap.extract_top();
            Swap bestOneMove = {nextBestVertex, NetworKit::none, improvementValues[nextBestVertex]};
            Swap bestSwap = {nextBestVertex, NetworKit::none, -std::numeric_limits<double>::max()};
            for (auto [neighbor, weight] : originalGraph_->weightNeighborRange(nextBestVertex)) {
                if (!heap.contains(neighbor))
                    continue;
                auto posGain = improvementValues[nextBestVertex] + improvementValues[neighbor]
                               + (partition.inA[nextBestVertex] ^ partition.inA[neighbor] ? 2.0 : -2.0) * weight;
                if (posGain > bestSwap.gain) {
                    bestSwap.v = neighbor;
                    bestSwap.gain = posGain;
                }
            }
            if ((bestSwap.v != NetworKit::none) && (bestSwap.gain > bestOneMove.gain)) {
                swapSuggestions.emplace_back(bestSwap);
                assert(heap.contains(bestSwap.v));
                heap.remove(bestSwap.v);
                updateImprovementValues(improvementValues, partition, bestSwap.u);
                updateImprovementValues(improvementValues, partition, bestSwap.v);
                heap.update_all();
            } else {
                swapSuggestions.emplace_back(bestOneMove);
                updateImprovementValues(improvementValues, partition, bestOneMove.u);
                heap.update_all();
            }

            numInnerIterations += 1;
            break;
        }

        auto [pos, maxGain] = maxGainSumIndex(swapSuggestions);
        if (gtEps(maxGain, 0.0, 10e-6)) {
#ifndef NDEBUG
            auto solValBefore = getPrimalSolutionValue(partition.inA);
            std::vector<bool> gotSwapped(originalGraph_->upperNodeIdBound());
#endif
            for (size_t i = 0; i <= pos; i++) {
                assert(partition.valid());
                auto &currentSwap = swapSuggestions[i];
                assert(currentSwap.u != NetworKit::none);
                if (currentSwap.v != NetworKit::none) {
                    if (partition.inA[currentSwap.u] != partition.inA[currentSwap.v]) {
                        partition.swapBetweenPartitions(currentSwap.u, currentSwap.v);
#ifndef NDEBUG
                        assert(!gotSwapped[currentSwap.u]);
                        assert(!gotSwapped[currentSwap.v]);
                        gotSwapped[currentSwap.u] = true;
                        gotSwapped[currentSwap.v] = true;
#endif
                    } else {
                        partition.swapToOtherPartition(currentSwap.u);
                        partition.swapToOtherPartition(currentSwap.v);
#ifndef NDEBUG
                        assert(!gotSwapped[currentSwap.u]);
                        gotSwapped[currentSwap.v] = true;
#endif
                    }
                } else {
                    partition.swapToOtherPartition(currentSwap.u);
                }
            }
#ifndef NDEBUG
            auto solValAfter = getPrimalSolutionValue(partition.inA);
            assert(equalAbsEps(solValBefore + maxGain, solValAfter, 10e-6));
#endif
            improvement = true;
        } else {
            improvement = false;
        }
        numOuterIterations += 1;
    }
    phase1Iters_ = numOuterIterations;
}

void KLHeuristic::phase2Optimization() {
    auto currentSolutionValue = getPrimalSolutionValue(solution_);

    std::vector<NetworKit::node> sideA;
    std::vector<NetworKit::node> sideB;
    for (auto node : originalGraph_->nodeRange()) {
        if (solution_[node])
            sideA.push_back(node);
        else
            sideB.push_back(node);
    }

    auto solutionA = phase1OnOneSide(sideA);
    auto solutionB = phase1OnOneSide(sideB);

    // if a node is fixed in A (or B) it will be fixed to A1 (or B1) in phase1OnOneSide,
    // so we only have to check the cases when A1 or B1 contains fixed nodes

    bool a1Fixed = false;
    bool b1Fixed = false;

    for (auto node : originalGraph_->nodeRange()) {
        assert(!(solution_[node] && solutionB[node]));  // node can not be in side A and B
        assert(!(!solution_[node] && solutionA[node])); // node can not be in side A and B
        if (vertexIsFixed_[node]) {
            assert(!(solution_[node] && !solutionA[node]));  // no fixed node should be in A2
            assert(!(!solution_[node] && !solutionB[node])); // no fixed node should be in B2
            if (solutionA[node])
                a1Fixed = true;
            else if (solutionB[node])
                b1Fixed = true;
        }
    }

    std::vector<bool> solutionA1B1 = solutionA; // sets A1 and B1 become A, A2 and B2 become B
    std::vector<bool> solutionA1B2 = solutionA; // sets A1 and B2 become A, A2 and B1 become B

    for (auto node : originalGraph_->nodeRange()) {
        if (solutionB[node]) // node is in B1
            solutionA1B1[node] = true;
        else if (!solution_[node]) // node is in B2
            solutionA1B2[node] = true;
    }

    auto solutionA1B1Value = getPrimalSolutionValue(solutionA1B1);
    auto solutionA1B2Value = getPrimalSolutionValue(solutionA1B2);

    if (!b1Fixed) { // B1 is not fixed, A1 may be fixed
        if (solutionA1B1Value >= currentSolutionValue && solutionA1B1Value >= solutionA1B2Value)
            solution_ = solutionA1B1;
        else if (solutionA1B2Value >= currentSolutionValue)
            solution_ = solutionA1B2;
    } else if (!a1Fixed) { // B1 is fixed, A1 is not fixed
        if (solutionA1B2Value >= currentSolutionValue && solutionA1B2Value >= solutionA1B1Value)
            solution_ = solutionA1B2;
        else if (solutionA1B1Value >= currentSolutionValue)
            solution_ = negateBoolVector(solutionA1B1);
    } else { // A1 and B1 are fixed
        if (solutionA1B2Value >= currentSolutionValue)
            solution_ = solutionA1B2;
    }
}

std::vector<bool> KLHeuristic::phase1OnOneSide(const std::vector<NetworKit::node> &nodes) const {
    std::vector<bool> solution(originalGraph_->numberOfNodes(), false);
    if (nodes.empty()) {
        return solution;
    } else if (nodes.size() == 1) {
        solution[nodes[0]] = true; // fix node to a1/b1 in case its fixed
        return solution;
    }
    auto subgraph = NetworKit::GraphTools::subgraphFromNodes(*originalGraph_, nodes.begin(), nodes.end(), false);
    auto nodeIdMap = NetworKit::GraphTools::getContinuousNodeIds(subgraph);
    auto compactSubgraph = NetworKit::GraphTools::getCompactedGraph(subgraph, nodeIdMap);

    std::vector<NetworKit::node> fixedVertices;
    std::vector<bool> initSolution(compactSubgraph.numberOfNodes(), false);
    for (auto node : nodes) {
        assert(subgraph.hasNode(node));
        if (vertexIsFixed_[node]) {
            assert(compactSubgraph.hasNode(nodeIdMap.at(node)));
            fixedVertices.push_back(nodeIdMap.at(node));
            initSolution[nodeIdMap.at(node)] = true;
        }
    }

    KLHeuristic klHeuristic(&compactSubgraph);
    klHeuristic.phase1Optimization(initSolution, fixedVertices);

    auto compactSolution = klHeuristic.getPrimalSolution();
    for (auto node : subgraph.nodeRange()) {
        if (compactSolution[nodeIdMap.at(node)]) {
            solution[node] = true;
        }
    }
    return solution;
}

std::vector<NetworKit::edgeweight> KLHeuristic::improvementValuesFromScratch(const BiPartition &partition) const {
    auto improvementValues = std::vector<NetworKit::edgeweight>(partition.inA.size(), 0.0);

    for (auto u : originalGraph_->nodeRange()) {
        bool sideU = partition.inA[u];
        for (auto [neighbor, weight] : originalGraph_->weightNeighborRange(u)) {
            bool neighborSide = partition.inA[neighbor];
            improvementValues[u] += (sideU ^ neighborSide ? -1.0 : 1.0) * weight;
        }
    }

    return improvementValues;
}

void KLHeuristic::updateImprovementValues(std::vector<NetworKit::edgeweight> &improvementValues, BiPartition &partition,
                                          NetworKit::node u) {
    assert(u < partition.inA.size());
    for (auto [neighbor, weight] : originalGraph_->weightNeighborRange(u)) {
        bool wasExternal = partition.inA[u] ^ partition.inA[neighbor];
        improvementValues[neighbor] += (wasExternal ? 2.0 : -2.0) * weight;
    }
}

std::tuple<size_t, NetworKit::edgeweight> KLHeuristic::maxGainSumIndex(const std::vector<Swap> &gainPairs) const {
    assert(!gainPairs.empty());
    NetworKit::edgeweight gainSum = gainPairs[0].gain;
    NetworKit::edgeweight maxGainSum = gainSum;
    size_t maxGainSumIndex = 0;

    for (size_t i = 1; i < gainPairs.size(); i++) {
        gainSum += gainPairs[i].gain;
        if (gainSum > maxGainSum) {
            maxGainSum = gainSum;
            maxGainSumIndex = i;
        }
    }

    return {maxGainSumIndex, maxGainSum};
}

bool KLHeuristic::belowMaxRuntime() const {
    return getTimeElapsed() < maxRuntime_;
}

bool KLHeuristic::belowDualBound(const std::vector<bool> &inA) const {
    return getPrimalSolutionValue(inA) < dualBound_;
}

std::chrono::milliseconds KLHeuristic::getTimeElapsed() const {
    return t_.elapsed();
}

double KLHeuristic::getPrimalSolutionValue(const std::vector<bool> &solution) const {
    return getCutValue(*originalGraph_, solution);
}

int KLHeuristic::getPhase1Iters() const {
    return phase1Iters_;
}
void KLHeuristic::fixVertices(const std::vector<NetworKit::node> &fixedVertices) {
    for (auto node : fixedVertices) {
        if (originalGraph_->hasNode(node))
            vertexIsFixed_[node] = true;
    }
}
std::vector<bool> KLHeuristic::negateBoolVector(std::vector<bool> &vector) const {
    for (auto node : originalGraph_->nodeRange()) {
        vector[node] = !vector[node];
    }
    return vector;
}

} // namespace sms