#include "sms/graph/graphs_statistics.hpp"

#include "sms/auxiliary/math.hpp"
#include "sms/auxiliary/statistics.hpp"
#include "sms/graph/bipartite.hpp"
#include "sms/graph/graphs.hpp"

namespace sms {

nlohmann::ordered_json graphStatsJson(const NetworKit::Graph &g, bool fast) {
    nlohmann::ordered_json j;

    j["num nodes"] = g.numberOfNodes();
    j["num edges"] = g.numberOfEdges();

    // Degrees
    std::vector<unsigned int> degrees;
    for (auto v : g.nodeRange()) {
        degrees.push_back(g.degree(v));
    }

    auto degreeStats = IntStatistics(degrees);
    j["degree stats"] = degreeStats.asJson();

    j["density"] =
        static_cast<double>(g.numberOfEdges()) * 2.0 / static_cast<double>(g.numberOfNodes() * (g.numberOfNodes() - 1));
    std::array<NetworKit::count, 4> numDegrees = {0, 0, 0, 0};
    for (auto u : g.nodeRange()) {
        if (g.degree(u) < 4)
            numDegrees[g.degree(u)] += 1;
    }
    j["num degree 0"] = numDegrees[0];
    j["num degree 1"] = numDegrees[1];
    j["num degree 2"] = numDegrees[2];
    j["num degree 3"] = numDegrees[3];

    // Weights
    bool isInt = true;
    unsigned int numZero = 0;
    std::vector<double> weights;
    double sumWeights = 0;
    double sumPositiveWeights = 0;
    double sumNegativeWeights = 0;
    for (auto e : g.edgeWeightRange()) {
        weights.push_back(e.weight);
        isInt &= isInteger(e.weight);
        if (e.weight == 0)
            numZero++;
        if (e.weight > 0)
            sumPositiveWeights += e.weight;
        else
            sumNegativeWeights += e.weight;
        sumWeights += e.weight;
    }

    auto weightStats = FloatStatistics(weights);
    j["weight stats"] = weightStats.asJson();
    if (isInt && (sumPositiveWeights - sumNegativeWeights) > maxSafeInteger<double>()) {
        j["weight stats"]["warning"] = "Sum of weights is too large for integer representation.";
        std::cerr << "All weights are integer, but numerical precision is not guaranteed." << std::endl;
    }

    if (isInt) {
        j["weight stats"]["gcd"] = edgeWeightDivisor(g);
    }

    j["weight type"] = isInt ? "integer" : "floating point";
    j["num zero weight edges"] = numZero;
    j["sum of weights"] = sumWeights;
    j["sum of positive weights"] = sumPositiveWeights;

    Bipartite bipartite(g);
    bipartite.run(false);
    j["bipartite"] = bipartite.isBipartiteGraph();
    bipartite.run(true);
    j["perfect cut"] = bipartite.isBipartiteGraph();

    // Components
    NetworKit::ConnectedComponents comps(g);
    comps.run();

    unsigned int numOneComp = 0;
    std::vector<unsigned int> componentSizes;

    for (const auto &c : comps.getComponents()) {
        unsigned int size = c.size();

        if (size > 1)
            componentSizes.push_back(size);
        else
            numOneComp++;
    }

    auto compStats = IntStatistics(componentSizes);

    j["number of cc"] = comps.getComponents().size();
    j["size of largest cc"] = compStats.max;
    j["size of smallest cc >1"] = compStats.min;
    j["mode of cc sizes"] = compStats.mode;
    j["number of size one cc"] = numOneComp;

    NetworKit::BiconnectedComponents biCons(g);
    biCons.run();

    unsigned int numTwoBiCons = 0;
    std::vector<unsigned int> biConsSizes;

    for (const auto &c : biCons.getComponents()) {
        unsigned int size = c.size();

        if (size > 2)
            biConsSizes.push_back(size);
        else if (size == 2)
            numTwoBiCons++;
    }

    j["number of bcc"] = biCons.getComponents().size();

    if (not biConsSizes.empty()) {
        auto biconStats = IntStatistics(biConsSizes);
        j["size of largest bcc"] = biconStats.max;
        j["size of smallest bcc >2"] = biconStats.min;
        j["mode of bcc sizes"] = biconStats.mode;
        j["num size two bcc"] = numTwoBiCons;
    }

    if (not fast) { // short chordless cycles
        SmallChordlessCycles ccs(g);
        ccs.run();
        j["num triangles"] = ccs.triangles.size();
        j["num 4-holes"] = ccs.fourHoles.size();

        auto mc = NetworKit::MaximalCliques(g);
        mc.run();
        auto cliques = mc.getCliques();

        std::vector<unsigned int> cliqueSizes;

        for (const auto &i : cliques) {
            if (i.size() >= 5)
                cliqueSizes.push_back(i.size());
        }

        if (!cliqueSizes.empty()) {
            j["size biggest clique found of size >= 5"] = *std::max_element(cliqueSizes.begin(), cliqueSizes.end());
        }
    }

    int stopIfLarger = 10;
    auto degen = degeneracy(g, stopIfLarger);
    if (degen.has_value()) {
        j["degeneracy"] = *degen;
    }

    return j;
}

void outputGraphStats(const NetworKit::Graph &g, std::ostream &out, bool fast) {
    out << "---------- Graph statistics -----------------" << std::endl;
    auto j = graphStatsJson(g, fast);
    out << j.dump(4) << std::endl;
    out << "---------------------------------------------" << std::endl;
}

std::optional<int> degeneracy(const NetworKit::Graph &g, unsigned int stopIfLarger) {
    if (g.numberOfNodes() < 1)
        throw std::runtime_error("Graph has no nodes");
    if (g.numberOfSelfLoops() > 0 || g.isDirected())
        throw std::runtime_error("Graph has self loops or is directed");

    std::vector<std::vector<NetworKit::node>> degreeBuckets(stopIfLarger + 1);
    std::vector<std::size_t> posInBucket(g.upperNodeIdBound(), 0);
    std::vector<NetworKit::count> degrees(g.upperNodeIdBound(), 0);

    NetworKit::count minDegree = std::numeric_limits<NetworKit::count>::max();
    for (auto u : g.nodeRange()) {
        NetworKit::count degree = g.degree(u);
        degrees[u] = degree;
        minDegree = std::min(minDegree, degree);
        if (degree <= stopIfLarger) {
            degreeBuckets[degree].push_back(u);
            posInBucket[u] = degreeBuckets[degree].size() - 1;
        }
    }

    if (minDegree == 0)
        return 0;

    std::vector<bool> active(g.upperNodeIdBound(), true);
    int numActiveNodes = g.numberOfNodes();
    NetworKit::count currentDegen = minDegree;
    while ((currentDegen <= stopIfLarger)) {
        if (minDegree == 0) {
            if (static_cast<std::size_t>(numActiveNodes) == degreeBuckets[0].size())
                break;
            for (unsigned int i = 1; i < degreeBuckets.size(); ++i) {
                if (!degreeBuckets[i].empty()) {
                    minDegree = i;
                    currentDegen = std::max(currentDegen, minDegree);
                    break;
                }
            }
            if (minDegree == 0)
                return std::nullopt;
        }
        if (degreeBuckets[minDegree].empty()) {
            minDegree++;
            currentDegen = std::max(currentDegen, minDegree);
            continue;
        }
        const auto nextNode = degreeBuckets[minDegree].back();
        degreeBuckets[minDegree].pop_back();
        active[nextNode] = false;
        numActiveNodes--;
        assert(degrees[nextNode] == minDegree);

        for (auto neighbor : g.neighborRange(nextNode)) {
            if (!active[neighbor])
                continue;
            const auto neighborDegree = degrees[neighbor];
            if (neighborDegree - 1
                > stopIfLarger) { // by reducing the degree by 1, the node still doesn't enter a bucket
                degrees[neighbor] = neighborDegree - 1;
                continue;
            } else if (neighborDegree <= stopIfLarger) { // the node is already in a bucket and needs to be removed
                auto pos = posInBucket[neighbor];
                // remove neighbor from bucket i
                auto lastNodeInBucket = degreeBuckets[neighborDegree].back();
                degreeBuckets[neighborDegree].pop_back();
                if (lastNodeInBucket != neighbor) {
                    degreeBuckets[neighborDegree][pos] = lastNodeInBucket;
                    posInBucket[lastNodeInBucket] = pos;
                }
            }

            // if degree was stopIfLarger, the node wasn't in a bucket before but gets added to one now
            // add neighbor to bucket i-1
            degreeBuckets[neighborDegree - 1].push_back(neighbor);
            degrees[neighbor] = neighborDegree - 1;
            posInBucket[neighbor] = degreeBuckets[neighborDegree - 1].size() - 1;
        }
        if (!degreeBuckets[minDegree - 1].empty()) { // minDegree can decrease by a maximum of 1
            minDegree--;
        } else if (degreeBuckets[minDegree].empty()) {
            minDegree++;
        }
        currentDegen = std::max(currentDegen, minDegree);
    }

    if (currentDegen > stopIfLarger)
        return {};
    else
        return currentDegen;
}

} // namespace sms
