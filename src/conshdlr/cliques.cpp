#include "sms/conshdlr/cliques.hpp"

#include <cmath>

#include "sms/auxiliary/bipartition.hpp"
#include "sms/auxiliary/math.hpp"
#include "sms/auxiliary/subset_enumerator.hpp"
#include "sms/instance/mc_solution.hpp"
#include "sms/pheur/kl_heuristic.hpp"
#include "sms/solver/enumeration_solver.hpp"

namespace sms {

struct weightedNode {
    NetworKit::node u;
    double weight;
};

CliqueSeparator::CliqueSeparator(const NetworKit::Graph &originalGraph) : lpWeightedGraph_(originalGraph) {}

void CliqueSeparator::updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
    // Guarantee numerical stability at 0 and 1
    NetworKit::edgeweight fixedVal = std::min(std::max(0., w), 1.);
    lpWeightedGraph_.setWeight(u, v, fixedVal);
}

const std::vector<NetworKit::node> &CliqueSeparator::getClique(unsigned int pos) const {
    return maximalCliques_[pos];
}

std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
CliqueSeparator::checkSimpleViolation(std::vector<NetworKit::node> clique) const {
    assert(clique.size() % 2 == 1);
    auto rhs = (clique.size() / 2) * (clique.size() + 1) / 2; // div 2 floor times ceil

    bool probablyNoViolation = true;

    NetworKit::Graph g(clique.size(), true, false);
    double offset = 0;
    for (unsigned int i = 0; i < clique.size(); i++) {
        for (unsigned int j = i + 1; j < clique.size(); j++) {
            auto u = clique[i];
            auto v = clique[j];
            assert(lpWeightedGraph_.hasEdge(u, v));
            NetworKit::edgeweight weight = lpWeightedGraph_.weight(u, v);
            probablyNoViolation &= equalEps(weight, 0.5, 0.1);
            g.addEdge(i, j, 1 - 2 * weight);
            offset += weight;
        }
    }

    if (probablyNoViolation)
        return {};

    if (clique.size() < enumerationLimit_) {
        EnumerationSolver es(g);
        es.run();

        if (es.bestSolutionValue() + offset > static_cast<double>(rhs)) {
            auto sol = es.getBestSolution();
            return bipartitionToIndexSets(sol, clique);
        }
    } else {
        KLHeuristic klHeuristic(&g);
        klHeuristic.phase1Optimization(std::vector<bool>(g.numberOfNodes(), false));
        auto sol = klHeuristic.getPrimalSolution();
        double solValue = solutionValue(g, sol);
        if (solValue + offset > static_cast<double>(rhs)) {
            return bipartitionToIndexSets(sol, clique);
        }
    }

    return {};
}

std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>>
CliqueSeparator::checkViolation(const std::vector<NetworKit::node> &clique, int lowerSubCliqueSize,
                                int upperSubCliqueSize, bool recursive) const {
    int cliqueSize = clique.size();

    if (cliqueSize < lowerSubCliqueSize)
        return {};
    if (static_cast<unsigned int>(cliqueSize) > smallCliqueSize_)
        return {checkBigClique(clique)};

    // all vertex subsets of the clique of size i with at least 5 vertices
    lowerSubCliqueSize = std::max(lowerSubCliqueSize, 5);
    upperSubCliqueSize = std::min(cliqueSize, upperSubCliqueSize);

    lowerSubCliqueSize = lowerSubCliqueSize % 2 == 0 ? lowerSubCliqueSize + 1 : lowerSubCliqueSize;
    upperSubCliqueSize = upperSubCliqueSize % 2 == 0 ? upperSubCliqueSize - 1 : upperSubCliqueSize;

    std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>> violations;

    if (fromSmallToBig_) {
        for (int i = lowerSubCliqueSize; i <= upperSubCliqueSize; i += 2) {
            violations = checkSubCliques(clique, i, lowerSubCliqueSize, upperSubCliqueSize, recursive);
            if (!violations.empty())
                return violations;
        }
    } else {
        for (int i = upperSubCliqueSize; i >= lowerSubCliqueSize; i -= 2) {
            violations = checkSubCliques(clique, i, lowerSubCliqueSize, upperSubCliqueSize, recursive);
            if (!violations.empty())
                return violations;
        }
    }
    return {};
}

// simple check of full clique
std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
CliqueSeparator::checkBigOddCliqueSimple(const std::vector<NetworKit::node> &clique) const {
    assert(clique.size() % 2 == 1);

    bool probablyNoViolation = false;

    double sum = 0.;
    NetworKit::Graph g(clique.size(), true, false);
    for (unsigned int i = 0; i < clique.size(); i++) {
        for (unsigned int j = i + 1; j < clique.size(); j++) {
            auto u = clique[i];
            auto v = clique[j];
            assert(lpWeightedGraph_.hasEdge(u, v));
            NetworKit::edgeweight weight = lpWeightedGraph_.weight(u, v);
            probablyNoViolation |= equalEps(weight, 0.5, 0.1);
            g.addEdge(i, j, 1 - 2 * weight);
            sum += weight;
        }
    }

    auto rhs = static_cast<double>((clique.size() / 2) * ((clique.size() + 1) / 2)); // div 2 floor times ceil
    if (sum > rhs) {
        return {clique, {}};
    }

    // only search for violations if chance of violation seems high
    if (probablyNoViolation)
        return {};

    KLHeuristic klHeuristic(&g);
    klHeuristic.phase1Optimization(std::vector<bool>(g.numberOfNodes(), false));

    auto sol = klHeuristic.getPrimalSolution();
    double solValue = 0;
    for (auto edge : g.edgeWeightRange()) {
        solValue += (sol[edge.u] ^ sol[edge.v]) * edge.weight;
        if (solValue + sum > rhs) {
            return bipartitionToIndexSets(sol, clique);
        }
    }
    return {};
}

// simple check of full clique
std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
CliqueSeparator::checkBigOddClique(const std::vector<NetworKit::node> &clique) const {
    assert(clique.size() % 2 == 1);

    auto res = checkBigOddCliqueSimple(clique);
    if (!res.first.empty() || !res.second.empty()) {
        return res;
    } else { // no violated constrain was found -> tripartition
        return checkTriPartition(clique);
    }
}

// simple check of full clique
std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
CliqueSeparator::checkBigClique(const std::vector<NetworKit::node> &clique) const {
    assert(clique.size() > 2);
    if (clique.size() % 2 == 1) {
        return checkBigOddClique(clique);
    }

    NetworKit::node missing = clique[0];
    auto subClique = std::vector<NetworKit::node>(clique.begin() + 1, clique.end());
    auto res = checkBigOddClique(subClique);
    if (!res.first.empty() || !res.second.empty()) {
        return res;
    }

    for (unsigned int i = 0; i < clique.size() - 1; i++) {
        std::swap(missing, subClique[i]);
        res = checkBigOddClique(subClique);
        if (!res.first.empty() || !res.second.empty()) {
            return res;
        }
    }

    return {};
}

void CliqueSeparator::sortCliquesBySize() {
    std::sort(maximalCliques_.begin(), maximalCliques_.end(),
              [](const auto &a, const auto &b) { return a.size() > b.size(); });
}

std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
CliqueSeparator::checkTriPartition(const std::vector<NetworKit::node> &clique) const {
    assert(clique.size() > 3);
    assert(clique.size() % 2 == 1);

    std::vector<weightedNode> shiftedWeightedDegrees(clique.size());

    for (unsigned i = 0; i < clique.size(); ++i) {
        shiftedWeightedDegrees[i] = {clique[i], 0.};
    }

    for (unsigned i = 0; i < clique.size(); i++) {
        for (unsigned j = i + 1; j < clique.size(); j++) {
            auto weight = lpWeightedGraph_.weight(clique[i], clique[j]);
            shiftedWeightedDegrees[i].weight += std::abs(0.5 - weight);
            shiftedWeightedDegrees[j].weight += std::abs(0.5 - weight);
        }
    }

    std::sort(shiftedWeightedDegrees.begin(), shiftedWeightedDegrees.end(),
              [](const weightedNode &a, const weightedNode &b) { return a.weight < b.weight; });

    std::vector<NetworKit::node> subClique = {};

    for (int i = 0; i < 5; i++) {
        subClique.push_back(shiftedWeightedDegrees.back().u);
        shiftedWeightedDegrees.pop_back();
    }

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>> res;

    while (subClique.size() < clique.size()) {
        if (subClique.size() <= smallCliqueSize_) {
            res = checkSimpleViolation(subClique);
        } else {
            res = checkBigOddCliqueSimple(subClique);
        }
        if (!res.first.empty() || !res.second.empty()) {
            return res;
        }

        // no violated constraint found -> increase subClique by 2 nodes
        for (int i = 0; i < 2; i++) {
            subClique.push_back(shiftedWeightedDegrees.back().u);
            shiftedWeightedDegrees.pop_back();
        }
    }
    return {};
}

std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>>
CliqueSeparator::checkSubCliques(const std::vector<NetworKit::node> &clique, int subCliqueSize, int lowerSubCliqueSize,
                                 int upperSubCliqueSize, bool recursive) const {
    if (subCliqueSize > upperSubCliqueSize || subCliqueSize < lowerSubCliqueSize)
        return {};

    int cliqueSize = clique.size();
    assert(cliqueSize >= subCliqueSize);
    assert(subCliqueSize >= 5 && subCliqueSize % 2 == 1);

    std::vector<std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>> violations;

    // all vertex subsets of the clique of size subCliqueSize
    auto se = sms::SubsetEnumerator(cliqueSize, subCliqueSize);
    while (!se.done()) {
        // build clique subset
        std::vector<NetworKit::node> subset;
        subset.reserve(subCliqueSize);
        std::vector<NetworKit::node> remainingClique;
        remainingClique.reserve(cliqueSize - subCliqueSize);
        for (int k = 0; k < cliqueSize; k++) {
            auto &appendTo = se.contains(k) ? subset : remainingClique;
            appendTo.push_back(clique[k]);
        }

        // call clique function
        auto res = checkSimpleViolation(subset);
        if (!res.first.empty() || !res.second.empty()) {
            violations.push_back(res);
            if (recursive) {
                lowerSubCliqueSize = fromSmallToBig_ ? subCliqueSize : lowerSubCliqueSize;
                upperSubCliqueSize = fromSmallToBig_ ? upperSubCliqueSize : subCliqueSize;
                if (std::ssize(remainingClique) >= lowerSubCliqueSize
                    && std::ssize(remainingClique) <= upperSubCliqueSize) {
                    auto violatedPairs =
                        checkViolation(remainingClique, lowerSubCliqueSize, upperSubCliqueSize, recursive);
                    for (const auto &violatedPair : violatedPairs) {
                        violations.push_back(violatedPair);
                    }
                }
            }
            return violations;
        }
        se.next();
    }
    return {};
}

} // namespace sms
