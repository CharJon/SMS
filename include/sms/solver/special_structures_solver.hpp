#ifndef SMS_SPECIAL_STRUCTURES_SOLVER_HPP
#define SMS_SPECIAL_STRUCTURES_SOLVER_HPP

#include "abstract_solver.hpp"
#include <optional>

#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace sms {

/**
 * Special graphs / graph structures for which MaxCut is solvable in polynomial time:
 * - Bipartite with only positive edges
 * - Graphs for which all positive edges form a bipartite graph and all negative edges are within one of the partitions
 * - Complete graphs where all edges have the same positive weight
 * - Graphs which contain a K_n/2_n/2 subgraphs with only positive edges of the same weight
 */
class SpecialStructuresSolver : public MaxCutSolver {
public:
    explicit SpecialStructuresSolver(const NetworKit::Graph &g);

    bool applicable() override;

    void run() override;

    std::vector<bool> getBestSolution() const override;

    bool optimalityProven() const override;

    void setTryFullKhalfhalfSearch(bool tryFullKhalfhalfSearch);

private:
    bool applicabilityChecked_ = false;
    bool applicable_ = false;

    // config
    bool tryFullKhalfhalfSearch_ = true;

    // Graph properties
    double positiveEdgesUnitWeight_ = -1;
    bool positiveEdgesEqualWeight_ = true;
    unsigned int numPositiveEdges_ = 0;
    unsigned int numNegativeEdges_ = 0;
    double sumPositiveEdgeWeights_ = 0;
    NetworKit::count n_ = 0, m_ = 0, numPossibleEdges_ = 0, numEdgesKHalfHalf_ = 0;

    // data
    std::vector<uint8_t> bestPartition_;

    // Computes properties of the graph needed to decide which solvers to run
    void computeGraphProperties(const NetworKit::Graph &g);

    // If there is a partition of the vertices such that all positive edges are in the cut
    // and all negative edges aren't then this finds the cut
    std::optional<std::vector<uint8_t>> getPositiveBipartitionNegativeSeparation(const NetworKit::Graph &g) const;

    // Computes a complement graph of G_+ = (V, E_+) where E_+ contains all positive weight edges of G
    NetworKit::Graph getPositiveEdgesOnlyComplementGraph(const NetworKit::Graph &g) const;

    // If the graph contains a K_{n/2,n/2} subgraph with only positive edges this finds the partition
    // (only useful if all positive edges have equal weight)
    std::optional<std::pair<std::vector<uint8_t>, double>>
    getPositiveKhalfhalfPartition(const NetworKit::Graph &g) const;

    // If all edges have negative weight this computes a max cut in poly time (via max flow)
    std::optional<std::pair<std::vector<uint8_t>, double>> getNegativeEdgesOnlyMaxCut(const NetworKit::Graph &g) const;

    // If only one edge has positive weight this computes a Max Cut
    std::optional<std::pair<std::vector<uint8_t>, double>>
    getOnlyOnePositiveEdgeMaxCut(const NetworKit::Graph &g) const;
};

} // namespace sms

#endif // SMS_SPECIAL_STRUCTURES_SOLVER_HPP
