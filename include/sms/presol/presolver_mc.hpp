#ifndef SMS_PRESOLVER_HPP
#define SMS_PRESOLVER_HPP

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"

#include "sms/auxiliary/bipartition.hpp"
#include "sms/graph/biconnected_partition.hpp"
#include "sms/graph/graphs.hpp"
#include "sms/presol/data_reducer_mc.hpp"
#include "sms/solver/special_structures_solver.hpp"

namespace sms {

/*
 * Decomposition and data reduction techniques
 */
class PresolverMC {

public:
    explicit PresolverMC(const NetworKit::Graph &graph, int emphasis = 1) : graph_(graph), emphasis_(emphasis) {
        if (graph_.numberOfNodes() != graph_.upperNodeIdBound())
            throw std::runtime_error("Graph is not compact");
        if (!graph_.isWeighted())
            throw std::runtime_error("Graph is not weighted");
        if (graph_.isDirected())
            throw std::runtime_error("Graph is directed");
        if (graph_.numberOfNodes() == 0)
            throw std::runtime_error("Graph is empty");
    }

    /*
     * @brief Run the presolver.
     * The result is a list of independent subgraphs
     */
    void run();

    unsigned int numberOfSubgraphs() const { return subgraphs_.size(); }

    unsigned int numberOfRemainingSubgraphs() const {
        int cnt = 0;
        for (const auto &r : reducers_) {
            const NetworKit::Graph &g = r.getGraph();
            if (g.numberOfNodes() > 0)
                cnt += 1;
        }
        return cnt;
    }

    nlohmann::ordered_json getStats() const;

    nlohmann::ordered_json getCompactStats() const;

    /*
     * @brief Add a solution for an independent reduced subgraph.
     * The solution has to assign every existing vertex to a partition
     */
    void addSolution(unsigned int subgraphPos, const std::vector<Partition> &sol) {
#ifndef NDEBUG
        const auto &redSubgraph = reducedSubgraph(subgraphPos);
        assert(sol.size() == redSubgraph.upperNodeIdBound());
        for (unsigned int i = 0; i < sol.size(); i++) {
            if (sol[i] == Partition::kUNASSIGNED)
                assert(!redSubgraph.hasNode(i));
            else
                assert(redSubgraph.hasNode(i));
        }
#endif
        if (subgraphPos < solutions_.size()) {
            solutions_[subgraphPos] = sol;
        } else {
            for (unsigned int i = 0; i < subgraphPos - solutions_.size() + 1; i++)
                solutions_.emplace_back();
            solutions_[subgraphPos] = sol;
        }
    }

    // Create a full solution from all subgraph solutions (which have to be added before)
    std::vector<bool> recoverSolution();

    void addWarmstartSolution(const std::vector<bool> &solution) {
        if (solution.size() != graph_.upperNodeIdBound())
            throw std::runtime_error("Invalid warmstart solution size");
        warmStartSolution_ = solution;
    }

    bool hasWarmstartSolution() const { return !warmStartSolution_.empty(); }

    // Return warmstart solution for a specific subgraph
    std::vector<bool> getSubgraphWarmstartSolution(int subgraphPos) const {
        if (warmStartSolution_.empty())
            throw std::runtime_error("No warmstart solution available");

        std::vector<bool> solution(reducedSubgraph(subgraphPos).upperNodeIdBound(), false);
        for (auto u : reducedSubgraph(subgraphPos).nodeRange()) {
            solution[u] = warmStartSolution_[subgraphVertexToOriginalVertex_[subgraphPos][u]];
        }

        return solution;
    }

    const std::vector<NetworKit::Graph> &subgraphs() const { return subgraphs_; }

    const NetworKit::Graph &reducedSubgraph(unsigned int i) const { return reducers_[i].getGraph(); }

    const std::vector<NetworKit::node> &roots() const { return roots_; }

    const std::vector<std::vector<NetworKit::node>> &toOrig() const { return subgraphVertexToOriginalVertex_; }

    std::chrono::milliseconds elapsed() const { return elapsed_; }

    int numRemainingVertices() const {
        if (!hasRun_)
            throw std::runtime_error("No call to run() yet!");
        int count = static_cast<int>(reducedSubgraph(0).numberOfNodes());
        for (unsigned int i = 1; i < subgraphs_.size(); i++) {
            // minus one, to not count articulation points twice
            // empty subgraphs are not counted
            count += std::max(static_cast<int>(reducedSubgraph(i).numberOfNodes()) - 1, 0);
        }
        return count;
    }

    int numRemainingEdges() const {
        if (!hasRun_)
            throw std::runtime_error("No call to run() yet!");
        int count = 0;
        for (unsigned int i = 0; i < subgraphs_.size(); i++) {
            count += static_cast<int>(reducedSubgraph(i).numberOfEdges());
        }
        return count;
    }

    void setOldPresolve();

    void setTriangleEmphasis(int newEmphasis) { triangleEmphasis_ = newEmphasis; }

private:
    const NetworKit::Graph &graph_;

    // The list of subgraphs is ordered (result of some DFS graph traversal)!
    std::vector<NetworKit::Graph> subgraphs_;
    // Root of corresponding subgraph
    std::vector<NetworKit::node> roots_; // root IDs in subgraph (not original graph)
    std::vector<std::vector<NetworKit::node>> subgraphVertexToOriginalVertex_;
    std::vector<std::vector<Partition>> solutions_;
    std::vector<bool> warmStartSolution_;
    std::vector<DataReducerMC> reducers_;
    std::chrono::milliseconds elapsed_;
    bool hasRun_ = false;
    const int emphasis_;
    bool oldPresolve_ = false;
    int triangleEmphasis_ = 3;

    void addComponent(const std::vector<NetworKit::node> &component, NetworKit::node root);

    void decompositionPhase();

    /*
     * Returns numer of components
     */
    unsigned int biconnectedComponentsDecomposition(const NetworKit::Graph &graph,
                                                    const std::vector<NetworKit::node> &toOrig);

    void dataReductionPhase();
};

} // namespace sms

#endif // SMS_PRESOLVER_HPP
