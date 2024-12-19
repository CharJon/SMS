#ifndef SMS_KL_HEURISTIC_HPP
#define SMS_KL_HEURISTIC_HPP

#include <cassert>
#include <chrono>
#include <utility>
#include <vector>

#include "networkit/auxiliary/Timer.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"

namespace sms {

struct Swap {
    NetworKit::node u;
    NetworKit::node v;
    NetworKit::edgeweight gain;
};

struct BiPartition {
    std::vector<NetworKit::node> A;
    std::vector<NetworKit::node> B;
    std::vector<bool> inA;
    std::vector<size_t> position;

    explicit BiPartition(const std::vector<bool> &sides) : inA(sides), position(inA.size(), 0) {
        for (size_t i = 0; i < inA.size(); i++) {
            if (inA[i]) {
                A.push_back(i);
                position[i] = A.size() - 1;
            } else {
                B.push_back(i);
                position[i] = B.size() - 1;
            }
        }
    }

    void swapBetweenPartitions(NetworKit::node u, NetworKit::node v) {
        assert(u != v);
        assert(inA[u] != inA[v]);
        auto &partitionOfU = inA[u] ? A : B;
        auto &partitionOfV = inA[v] ? A : B;
        std::swap(partitionOfU[position[u]], partitionOfV[position[v]]);
        std::swap(position[u], position[v]);
        std::vector<bool>::swap(inA[u], inA[v]);
    }

    void swapToOtherPartition(NetworKit::node u) {
        auto &partitionOfU = inA[u] ? A : B;
        auto &otherPartition = inA[u] ? B : A;
        partitionOfU[position[u]] = partitionOfU.back();
        position[partitionOfU.back()] = position[u];
        partitionOfU.pop_back();
        otherPartition.push_back(u);
        position[u] = otherPartition.size() - 1;
        inA[u] = !inA[u];
    }

    bool valid() {
        bool isValid = true;
        for (size_t i = 0; i < inA.size(); i++) {
            auto &partition = inA[i] ? A : B;
            isValid &= partition[position[i]] == i;
            assert(isValid);
        }
        return isValid;
    }
};

class KLHeuristic {
public:
    NetworKit::Graph const *const originalGraph_;

    explicit KLHeuristic(NetworKit::Graph const *const originalGraph,
                         double dualBound = std::numeric_limits<double>::max())
        : originalGraph_(originalGraph),
          solution_(originalGraph->upperNodeIdBound(), false),
          vertexIsFixed_(originalGraph->upperNodeIdBound(), false),
          dualBound_(dualBound) {
        if (originalGraph_->numberOfNodes() != originalGraph_->upperNodeIdBound()) {
            throw std::runtime_error("Graph needs to be compact!");
        }
    }

    void phase1Optimization(const std::vector<bool> &initSol, const std::vector<NetworKit::node> &fixedVertices = {});
    // optimise on a given solution with fixed vertices

    void phase2Optimization(); // optimize without solution and reoptimize anything found

    const std::vector<bool> &getPrimalSolution() const;

    std::vector<NetworKit::edgeweight> improvementValuesFromScratch(const BiPartition &partition) const;

private:
    // data
    std::vector<bool> solution_;
    std::vector<bool> vertexIsFixed_;
    Aux::Timer t_;

    // config
    std::chrono::seconds maxRuntime_ = std::chrono::seconds(3600);
    double dualBound_ = std::numeric_limits<double>::max();

    // stats
    int phase1Iters_ = 0;

    // start from given partition
    void phase1OptimizationSubRoutine(BiPartition &partition);

    // optimize one side
    std::vector<bool> phase1OnOneSide(const std::vector<NetworKit::node> &nodes) const;

    void updateImprovementValues(std::vector<NetworKit::edgeweight> &improvementValues, BiPartition &partition,
                                 NetworKit::node u);

    // Given a list of suggested swaps, find the end of the maximum partial sum
    std::tuple<size_t, NetworKit::edgeweight> maxGainSumIndex(const std::vector<Swap> &gainPairs) const;

    bool belowMaxRuntime() const;
    bool belowDualBound(const std::vector<bool> &) const;

    std::chrono::milliseconds getTimeElapsed() const;

    double getPrimalSolutionValue(const std::vector<bool> &) const;

    int getPhase1Iters() const;

    void fixVertices(const std::vector<NetworKit::node> &);

    std::vector<bool> negateBoolVector(std::vector<bool> &) const;
};

} // namespace sms

#endif // SMS_KL_HEURISTIC_HPP
