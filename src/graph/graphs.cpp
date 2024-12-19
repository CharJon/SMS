#include "sms/graph/graphs.hpp"

#include <map>
#include <vector>

#include "networkit/graph/Graph.hpp"

namespace sms {

NetworKit::Graph unweightedComplementGraph(const NetworKit::Graph &graph) {
    if (graph.isDirected() || (graph.numberOfNodes() != graph.upperNodeIdBound())) {
        throw std::runtime_error("Not implemented for directed graphs or non-compact graphs.");
    }

    auto complementGraph = NetworKit::Graph(graph.upperNodeIdBound(), false, false);

    auto lastNeighbor = std::vector<NetworKit::node>(
        graph.upperNodeIdBound(),
        std::numeric_limits<NetworKit::node>::max()); // max of unsigned type here allows us to get 0 if we add 1
    for (NetworKit::node u = 0; u < graph.upperNodeIdBound(); u++) {
        for (auto neighbor : graph.neighborRange(u)) {
            for (auto v = lastNeighbor[neighbor] + 1; v < std::min(u, neighbor); v++) {
                assert(!complementGraph.hasEdge(neighbor, v));
                complementGraph.addEdge(neighbor, v);
            }
            lastNeighbor[neighbor] = u;
        }
    }

    for (NetworKit::node u = 1; u < graph.upperNodeIdBound(); u++) {
        for (NetworKit::node v = lastNeighbor[u] + 1; v < u; v++) {
            complementGraph.addEdge(u, v);
        }
    }

    return complementGraph;
}

std::tuple<NetworKit::node, NetworKit::count> maxDegreeNode(const NetworKit::Graph &g) {
    assert(g.numberOfNodes() > 0);
    NetworKit::node maxDegreeNode = NetworKit::none;
    NetworKit::count maxDegree = 0;
    for (auto u : g.nodeRange()) {
        auto currentDegree = g.degree(u);
        if (currentDegree > maxDegree) {
            maxDegree = currentDegree;
            maxDegreeNode = u;
        }
    }
    assert(maxDegreeNode != NetworKit::none);
    return {maxDegreeNode, maxDegree};
}

CompactGraph compactGraph(const NetworKit::Graph &graph, int seed) {
    auto compactToOrig = std::vector<NetworKit::node>(graph.numberOfNodes());
    std::unordered_map<NetworKit::node, NetworKit::node> origToCompact;

    NetworKit::Graph compactGraph = NetworKit::Graph(graph.numberOfNodes(), graph.isWeighted(), graph.isDirected());
    NetworKit::node c = 0;

    if (seed > 0) {
        auto nodeIt = graph.nodeRange();
        std::vector<NetworKit::node> allVertices(nodeIt.begin(), nodeIt.end());
        std::mt19937 randomGen(seed);
        std::shuffle(allVertices.begin(), allVertices.end(), randomGen);
        for (auto v : allVertices) {
            compactToOrig[c] = v;
            origToCompact[v] = c;
            c++;
        }
    } else {
        for (auto v : graph.nodeRange()) {
            compactToOrig[c] = v;
            origToCompact[v] = c;
            c++;
        }
    }

    for (auto e : graph.edgeWeightRange()) {
        compactGraph.addEdge(origToCompact[e.u], origToCompact[e.v], e.weight);
    }

    return {std::move(compactGraph), std::move(compactToOrig), std::move(origToCompact)};
}

CompactGraph inducedSubgraphCompact(const NetworKit::Graph &g, const std::vector<NetworKit::node> &nodes) {
    std::unordered_map<NetworKit::node, NetworKit::node> originalToTransformed;
    uint64_t currentId = 0;
    for (auto u : nodes) {
        originalToTransformed[u] = currentId;
        currentId += 1;
    }

    auto inducedGraph = NetworKit::Graph(nodes.size(), g.isWeighted(), g.isDirected());

    if (g.isWeighted()) {
        for (auto v : nodes) {
            for (auto [neighbor, weight] : g.weightNeighborRange(v)) {
                if (g.isDirected() || ((v < neighbor) && originalToTransformed.count(neighbor))) {
                    inducedGraph.addEdge(originalToTransformed[v], originalToTransformed[neighbor], weight);
                }
            }
        }
    } else {
        for (auto v : nodes) {
            for (auto neighbor : g.neighborRange(v)) {
                if (g.isDirected() || ((v < neighbor) && originalToTransformed.count(neighbor))) {
                    inducedGraph.addEdge(originalToTransformed[v], originalToTransformed[neighbor]);
                }
            }
        }
    }

    return {std::move(inducedGraph), nodes, std::move(originalToTransformed)};
}

NetworKit::edgeweight edgeWeightDivisor(const NetworKit::Graph &graph) {
    if (graph.numberOfEdges() < 1)
        return 1.;

    bool intWeights = true;
    bool allAbsTheSame = true;
    NetworKit::edgeweight firstWeight = (*(graph.edgeWeightRange().begin())).weight;
    NetworKit::edgeweight minEdgeWeight = firstWeight;
    NetworKit::edgeweight maxEdgeWeight = firstWeight;
    auto gcd = static_cast<int64_t>(firstWeight);

    for (auto e : graph.edgeWeightRange()) {
        allAbsTheSame &= (std::fabs(e.weight) == std::fabs(firstWeight)) || (e.weight == 0);
        minEdgeWeight = std::min(minEdgeWeight, e.weight);
        maxEdgeWeight = std::max(maxEdgeWeight, e.weight);
        gcd = std::gcd(gcd, static_cast<int64_t>(e.weight));
        intWeights &= isInteger(e.weight);
    }

    if (allAbsTheSame || (minEdgeWeight == maxEdgeWeight)) {
        return std::fabs(maxEdgeWeight);
    } else if (intWeights) {
        return static_cast<NetworKit::edgeweight>(gcd);
    }
    return 1.;
}

NetworKit::edgeweight degreeBasedScaling(const NetworKit::Graph &graph) {
    auto nodeRange = graph.nodeRange();
    if (graph.numberOfEdges() < 1)
        return 1.;

    // collect degrees of vertices induced by all edges with odd weight
    std::vector<int> oddDegrees(graph.upperNodeIdBound(), 0);

    for (auto e : graph.edgeWeightRange()) {
        if (!isInteger(e.weight)) {
            return 1.;
        }
        int isOdd = static_cast<int64_t>(std::fabs(e.weight)) % 2;
        oddDegrees[e.u] += isOdd;
        oddDegrees[e.v] += isOdd;
    }

    if (std::all_of(nodeRange.begin(), nodeRange.end(), [&oddDegrees](auto u) { return oddDegrees[u] % 2 == 0; }))
        return 2.;

    return 1.;
}

NetworKit::edgeweight neighborhoodAlphaMarksSet(const NetworKit::Graph &g, NetworKit::node u, NetworKit::node v,
                                                bool ignoreEachOther, std::vector<NetworKit::edgeweight> &marks) {
    assert(!g.isDirected());
    assert(g.hasNode(u));
    assert(g.hasNode(v));
    assert(marks.size() == g.upperNodeIdBound());

    if (g.degree(u) != g.degree(v)) {
        return 0.0;
    }

    if (g.degree(u) == 1) {
        auto firstNeighborU = g.getIthNeighbor(u, 0);
        if (ignoreEachOther && (firstNeighborU == v))
            return 1.0;
        auto firstNeighborV = g.getIthNeighbor(v, 0);
        if (firstNeighborU != firstNeighborV)
            return 0.0;
        return g.getIthNeighborWeight(u, 0) / g.getIthNeighborWeight(v, 0);
    }

    NetworKit::edgeweight alpha = 0.0;
    for (unsigned int i = 0; i < g.degree(v); i++) {
        auto neighbor = g.getIthNeighbor(v, i);
        auto weight = g.getIthNeighborWeight(v, i);
        double marksWeight = marks[neighbor];
        assert(weight != 0);
        if (neighbor == u && ignoreEachOther)
            continue;
        if (neighbor == u || marksWeight == 0) {
            return 0.;
        }
        if (alpha == 0) // first alpha
            alpha = marksWeight / weight;
        if (marksWeight != alpha * weight) {
            return 0.;
        }
    }
    return alpha;
}

NetworKit::edgeweight neighborhoodAlpha(const NetworKit::Graph &g, NetworKit::node u, NetworKit::node v,
                                        bool ignoreEachOther, std::vector<NetworKit::edgeweight> &marks) {
    assert(!g.isDirected());
    assert(g.hasNode(u));
    assert(g.hasNode(v));
    assert(marks.size() == g.upperNodeIdBound());
    assert(std::all_of(marks.begin(), marks.end(), [](auto &m) { return m == 0; }));

    if (g.degree(u) != g.degree(v)) {
        return 0.0;
    }
    for (auto [neighbor, weight] : g.weightNeighborRange(u)) {
        marks[neighbor] = weight;
    }
    auto alpha = neighborhoodAlphaMarksSet(g, u, v, ignoreEachOther, marks);
    for (auto [neighbor, weight] : g.weightNeighborRange(u)) {
        marks[neighbor] = 0.0;
    }
    return alpha;
}

NetworKit::edgeweight getCutValue(const NetworKit::Graph &g, const std::vector<bool> &solution) {
    NetworKit::edgeweight value = 0;
    for (auto edge : g.edgeWeightRange()) {
        value += (solution[edge.u] ^ solution[edge.v]) * edge.weight;
    }
    return value;
}

bool integerWeightedDegree(const NetworKit::Graph &g) {
    for (auto node : g.nodeRange()) {
        auto weightedDegree = g.weightedDegree(node);
        if (!isInteger(weightedDegree)) {
            return false;
        }
    }
    return true;
}

bool weightSumSafe(const NetworKit::Graph &g, int sign) {
    double sum = 0.;
    for (auto edge : g.edgeWeightRange()) {
        auto absWeight = sign * edge.weight;
        if (absWeight >= 0)
            sum += absWeight;
    }
    return (sum <= sms::maxSafeInteger<double>());
}

double integerScalar(const NetworKit::Graph &g, double maxScalar, double eps) {
    if (g.numberOfEdges() == 0)
        return 1.0;

    auto absMaxEdge = *std::max_element(
        g.edgeWeightRange().begin(), g.edgeWeightRange().end(),
        [](NetworKit::WeightedEdge e, NetworKit::WeightedEdge f) { return std::abs(e.weight) < std::abs(f.weight); });
    double absMaxWeight = std::abs(absMaxEdge.weight);

    double s = maxSafeInteger<double>() / absMaxWeight;
    maxScalar = std::min(maxScalar, s);

    if (maxScalar < 1)
        return 0;

    // simple scalars
    const std::array<double, 9> simpleScalars{3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};

    // get minimal absolute non-zero weight
    double minWeight = std::numeric_limits<double>::max();
    for (auto e : g.edgeWeightRange()) {
        auto weight = e.weight;
        if (weight != 0)
            minWeight = std::min(minWeight, std::abs(weight));
    }

    if (minWeight == std::numeric_limits<double>::max()) // all weights are 0
        return 1.0;

    double bestScalar = std::numeric_limits<double>::max();
    double scalar;

    // test reciprocal of the smallest value and simple scalars as scalar
    for (int i = 0; i < 2; ++i) {
        scalar = ((i == 0) ? (1.0 / minWeight) : 1.0);
        for (auto e : g.edgeWeightRange()) {
            auto weight = std::abs(e.weight);
            if (weight == 0)
                continue;

            while (scalar <= maxScalar && (weight * scalar < 0.5 || !isInteger(weight * scalar))) {
                bool simpleScalarWorked = false;
                for (auto simpleScalar : simpleScalars) {
                    if (isInteger(weight * scalar * simpleScalar)) {
                        scalar *= simpleScalar;
                        simpleScalarWorked = true;
                        break;
                    }
                }
                if (!simpleScalarWorked)
                    scalar *= 2.0;
            }
        }
        if (scalar <= maxScalar) {
            bestScalar = std::min(scalar, bestScalar);
            if (i == 0 && equalEps(scalar, 1. / minWeight, eps))
                return bestScalar;
        }
    }

    /* convert each value into a rational number, calculate the greatest common divisor of the numerators
     * and the smallest common multiple of the denominators
     */
    long gcd = 1;
    long scm = 1;

    bool initialised = false;

    for (auto e : g.edgeWeightRange()) {
        if (e.weight == 0.0)
            continue;
        int maxDenom = static_cast<int>(maxScalar);
        if (maxScalar > std::numeric_limits<int>::max())
            maxDenom = std::numeric_limits<int>::max();
        auto [num, denom] = rationalRepresentation(e.weight, maxDenom);
        if (denom == -1) // no rational representation found
            return 0;
        if (!initialised && !equalEps(num, 0, eps) && !equalEps(denom, 0, eps)) {
            gcd = std::abs(num);
            scm = denom;
            initialised = true;
        } else if (!equalEps(num, 0, eps)) {
            assert(denom > 0);
            gcd = std::gcd(gcd, std::abs(num));
            scm *= denom / std::gcd(scm, denom);
        }
    }
    scalar = static_cast<double>(scm) / static_cast<double>(gcd);
    if (scalar < bestScalar)
        bestScalar = scalar;
    if (bestScalar <= maxScalar)
        return bestScalar;
    return 0;
}

NetworKit::Graph scaleToInt(NetworKit::Graph g, double scalar) {
    if (scalar == 0)
        return {};
    else if (scalar == 1)
        return g;
    double weight;
    for (auto e : g.edgeWeightRange()) {
        weight = e.weight * scalar;
        assert(sms::isNearlyInteger(weight, 10e-6 * scalar));
        g.setWeight(e.u, e.v, std::round(weight));
    }
    return g;
}

} // namespace sms
