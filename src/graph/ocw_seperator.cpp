#include "sms/graph/ocw_seperator.hpp"

#include <iomanip>

#include "sms/auxiliary/math.hpp"

namespace sms {

OcwSeparator::OcwSeparator(const NetworKit::Graph &originalGraph)
    : n_(originalGraph.numberOfNodes()),
      originalGraph_(originalGraph),
      weights_(originalGraph.upperEdgeIdBound(), 0.),
      currentSource_(NetworKit::none),
      blacklist_(originalGraph.numberOfNodes(), false),
      intermediateVertices_() {
    if (originalGraph.numberOfNodes() != originalGraph.upperNodeIdBound())
        throw std::runtime_error("Graph needs to be compact.");
    if (!originalGraph.hasEdgeIds() || originalGraph.isDirected())
        throw std::runtime_error("Passed graph has no edge ids / is directed.");
    distances_ = std::vector<NetworKit::edgeweight>(2 * n_ + 1, kInfDist_);
    hopDistances_ = std::vector<int>(2 * n_ + 1, 2 * n_ + 2);
    previous_ = std::vector<NetworKit::node>(2 * n_ + 1, NetworKit::none);
}

void OcwSeparator::updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
    assert(u < n_);
    assert(v < n_);
    NetworKit::edgeweight fixedVal = std::min(std::max(0., w), 1.);
    auto idx = originalGraph_.edgeId(u, v);
    weights_[idx] = fixedVal;
}

void OcwSeparator::updateEdgeWeight(NetworKit::index idx, NetworKit::edgeweight ew) {
    weights_[idx] = std::min(std::max(0., ew), 1.);
}

void OcwSeparator::runDijsktra() {
    NetworKit::node target = twinId(currentSource_);

    // reset data
    std::fill(distances_.begin(), distances_.end(), kInfDist_); // timestamping, which run created dist
    std::fill(hopDistances_.begin(), hopDistances_.end(), n_);  // timestamping, which run created dist
    // previous_ does not need to be reset
    intermediateVertices_.clear();

    auto distSmaller = [&distances = distances_, &hopDistances = hopDistances_](NetworKit::node u, NetworKit::node v) {
        return distances[u] == distances[v] ? hopDistances[u] < hopDistances[v] : distances[u] < distances[v];
    };
    tlx::d_ary_addressable_int_heap<NetworKit::node, 2, decltype(distSmaller)> heap(distSmaller);

    NetworKit::edgeweight bestDist = 1.0 - 1e-6;
    NetworKit::edgeweight bestHopDist = 0;

    // start
    distances_[currentSource_] = 0.;
    hopDistances_[currentSource_] = 0;
    heap.push(currentSource_);

    auto updatePaths = [&](NetworKit::node target, NetworKit::edgeweight newDist, NetworKit::edgeweight newHopDist,
                           NetworKit::node predecessor) {
        // update dist, if better
        if (distances_[target] == newDist ? newHopDist < hopDistances_[target] : newDist < distances_[target]) {
            distances_[target] = newDist;
            hopDistances_[target] = newHopDist;
            previous_[target] = predecessor;
            if (heap.contains(target))
                heap.update(target);
            else
                heap.push(target);

            auto targetTwin = twinId(target);
            auto distSourceNeighborTarget = distances_[target] + distances_[targetTwin];
            auto hopDistSourceNeighborTarget = hopDistances_[target] + hopDistances_[targetTwin];
            if (distSourceNeighborTarget < bestDist
                || ((distSourceNeighborTarget == bestDist) && (hopDistSourceNeighborTarget < bestHopDist))) {
                bestDist = distSourceNeighborTarget;
                bestHopDist = hopDistSourceNeighborTarget;
                intermediateVertices_.push_back(target);
            } else if ((distSourceNeighborTarget == bestDist) && (hopDistSourceNeighborTarget < 2 * bestHopDist)) {
                // not strictly better, but maybe also good
                intermediateVertices_.push_back(target);
            }
        }
    };

    int numberOfEdgesSeen = 0;
    // the dijkstra is implicitly bi-directional
    while ((!heap.empty())                                // work left
           && (heap.top() != target)                      // target not found (just to be sure)
           && ((2 * distances_[heap.top()] <= bestDist))) // same or better distance possible
    //               || ((2 * distances_[heap.top()] == bestDist) && (2 * hopDistances_[heap.top()] < bestHopDist)))
    {
        NetworKit::node current = heap.extract_top();
        NetworKit::node currentOriginal = originalId(current);

        for (unsigned int i = 0; i < originalGraph_.degree(currentOriginal); i++) {
            auto [neighbor, edgeId] = originalGraph_.getIthNeighborWithId(currentOriginal, i);
            if (blacklist_[neighbor])
                continue;
            numberOfEdgesSeen += 1;
            const double newHopDist = hopDistances_[current] + 1;

            const auto stayWeight = weights_[edgeId];
            const auto totalStayDist = distances_[current] + stayWeight;
            const auto stayNeighbor = toSide(neighbor, side(current));
            updatePaths(stayNeighbor, totalStayDist, newHopDist, current);
            const auto crossWeight = 1.0 - stayWeight;
            const auto totalCrossDist = distances_[current] + crossWeight;
            const auto crossNeighbor = toSide(neighbor, !side(current));
            updatePaths(crossNeighbor, totalCrossDist, newHopDist, current);
        }
    }
    if (intermediateVertices_.empty())
        blacklistIntReachable(originalId(currentSource_));
    std::sort(intermediateVertices_.begin(), intermediateVertices_.end(), [this](auto u, auto v) {
        return distances_[u] == distances_[v] ? hopDistances_[u] < hopDistances_[v] : distances_[u] < distances_[v];
    });
}

std::vector<ParityEdge> OcwSeparator::getShortestPath(const unsigned long bestIntermediate) {
    std::vector<ParityEdge> path;
    auto currentNode = bestIntermediate;
    path.reserve(hopDistances_[currentNode]);
    while (currentNode != currentSource_) {
        NetworKit::node prevNode = previous_[currentNode];
        NetworKit::node origPrevNode = originalId(prevNode);
        NetworKit::node origCurrentNode = originalId(currentNode);
        if (sameSide(currentNode, prevNode)) {
            path.emplace_back(origCurrentNode, origPrevNode, true);
        } else {
            path.emplace_back(origCurrentNode, origPrevNode, false);
        }
        currentNode = prevNode;
    }
    return path;
}

OddClosedWalk OcwSeparator::intermediateVertexToOcw(NetworKit::node intermediate) {
    auto path = getShortestPath(intermediate);
    assert(path.size() == static_cast<std::size_t>(hopDistances_[intermediate]));
    assert(path.size() != 0);

    auto twinPath = getShortestPath(twinId(intermediate));
    assert(twinPath.size() == static_cast<std::size_t>(hopDistances_[twinId(intermediate)]));

    if (twinPath.size() == 0) {
        assert(path.front().u == path.back().v);
        OddClosedWalk ocw(originalId(intermediate));
        for (unsigned int i = 0; i < path.size(); i++) {
            auto [start, end, stay] = path[i];
            ocw.addEdge(end, stay);
        }
        return ocw;
    }

    assert(twinPath.size() != 0);
    assert(path.front().u == twinPath.front().u);
    assert(path.back().v == twinPath.back().v);

    // path and twin path implicitly start at the same vertex (bestIntermediate)
    // find first position from back, where they differ
    std::size_t endPath = path.size() - 1;
    std::size_t endTwinPath = twinPath.size() - 1;
    while ((path[endPath].v == twinPath[endTwinPath].v) && (endPath > 0) && (endTwinPath > 0)) {
        endPath--;
        endTwinPath--;
        auto sharedVertex = path[endPath].v;
        blacklist_[sharedVertex] = true;
    }
    if (path[endPath].v != twinPath[endTwinPath].v) {
        endPath++;
        endTwinPath++;
    }
    assert(path[endPath].v == twinPath[endTwinPath].v);

    OddClosedWalk ocw(originalId(intermediate));
    // skip last, because it will be first for the twin path
    for (unsigned int i = 0; i <= endPath; i++) {
        auto [start, end, stay] = path[i];
        ocw.addEdge(end, stay);
    }
    for (int i = endTwinPath; i >= 0; i--) {
        auto [start, end, stay] = twinPath[i];
        ocw.addEdge(start, stay);
    }
    assert(ocw.isValid(originalGraph_));
    return ocw;
}

std::optional<OddClosedWalk> OcwSeparator::getMostViolatedOcw(NetworKit::node source) {
    currentSource_ = source;
    if (blacklist_[currentSource_])
        return {};
    runDijsktra();
    if (intermediateVertices_.empty()) {
        return {};
    }

    return intermediateVertexToOcw(intermediateVertices_.back());
}

std::vector<OddClosedWalk> OcwSeparator::simpleOddCycles(const OddClosedWalk &ocw) const {
    assert(ocw.isValid());
    std::vector<OddClosedWalk> simpleOddCycles = {};

    auto previousPosition = std::vector<int64_t>(n_, -1);
    int64_t lastLeftPos = -1;
    for (uint64_t i = 0; i < ocw.size() + 1; i++) {
        const auto cycleVertex = ocw.getIthNode(i);
        if (lastLeftPos < previousPosition[cycleVertex]) {
            lastLeftPos = previousPosition[cycleVertex];
            if (i - lastLeftPos > 2) {
                OddClosedWalk simpleOddClosedWalk = ocw.extract(lastLeftPos, i);
                assert(simpleOddClosedWalk.isSimple());
                assert(simpleOddClosedWalk.isValid());
                simpleOddCycles.push_back(simpleOddClosedWalk);
            }
        }
        previousPosition[cycleVertex] = i;
    }
    // if no cycle was found, the original cycle was simple already
    if (simpleOddCycles.empty()) {
        assert(ocw.isSimple());
        simpleOddCycles.push_back(ocw);
    }
    return simpleOddCycles;
}

std::vector<OddClosedWalk> OcwSeparator::extractCycles(const OddClosedWalk &ocw) const {
    assert(ocw.isValid());
    assert(ocw.isSimple());
    std::vector<OddClosedWalk> result = {};

    const uint64_t cycleLength = ocw.size();

    // precompute values
    auto numCrossingUpTo = std::vector<int>(cycleLength + 1, 0);
    auto pathLengthUpTo = std::vector<NetworKit::edgeweight>(cycleLength + 1, 0.);

    int numCrossings = 0;
    NetworKit::edgeweight pathLength = 0.;
    for (unsigned int i = 1; i <= cycleLength; i++) {
        auto curNode = ocw.getIthNode(i);
        auto prevNode = ocw.getIthNode(i - 1);
        numCrossings += (ocw.ithEdgeIsStay(i - 1) ? 0 : 1);
        auto edgeId = originalGraph_.edgeId(curNode, prevNode);
        pathLength += ocw.ithEdgeIsStay(i - 1) ? weights_[edgeId] : 1 - weights_[edgeId];
        numCrossingUpTo[i] = numCrossings;
        pathLengthUpTo[i] = pathLength;
    }

    // walk along cycle and find chords
    struct Chord {
        uint64_t start;
        uint64_t end;
        uint64_t length;
        bool crossing;
        bool inner;
    };
    std::vector<Chord> goodChords = {};

    for (uint64_t i = 2; i < cycleLength; i++) {
        auto curNode = ocw.getIthNode(i);
        for (auto cNeighbor : originalGraph_.neighborRange(curNode)) {
            // skip edges that already are part of cycle
            if (cNeighbor == ocw.getIthNode(i - 1) || (cNeighbor == ocw.getIthNode(i + 1)))
                continue;
            // find neighbor position in cycle, if edge forms a chord
            // MightDo (JC): This can be speed up by sorting nodes
            auto backwardsChordNeighborPos = NetworKit::none;
            for (uint64_t j = 0; j < i; j++) {
                if (ocw.getIthNode(j) == cNeighbor) {
                    backwardsChordNeighborPos = j;
                    break;
                }
            }
            if (backwardsChordNeighborPos != NetworKit::none) {
                auto j = backwardsChordNeighborPos;
                // check both possible cycles
                uint64_t innerCrossingDif = numCrossingUpTo[i] - numCrossingUpTo[j];
                double innerPathLength = pathLengthUpTo[i] - pathLengthUpTo[j];
                double outerPathLength = pathLengthUpTo[cycleLength] - innerPathLength;
                auto edgeId = originalGraph_.edgeId(ocw.getIthNode(i), ocw.getIthNode(j));
                NetworKit::edgeweight edgeLpValue = weights_[edgeId];

                bool innerIsEven = innerCrossingDif % 2 == 0;
                if (innerIsEven) {
                    if (innerPathLength - edgeLpValue < 0) {
                        goodChords.push_back({j, i, i - j + 1, true, true});
                    }
                    // outer is odd here
                    if (outerPathLength + edgeLpValue < 1) {
                        goodChords.push_back({j, i, cycleLength - i + j + 1, false, false});
                    }
                } else {
                    if (innerPathLength + edgeLpValue < 1) {
                        goodChords.push_back({j, i, i - j + 1, false, true});
                    }
                    // outer is even here
                    if (outerPathLength - edgeLpValue < 0) {
                        goodChords.push_back({j, i, cycleLength - i + j + 1, true, false});
                    }
                }
            }
        }
    }

    // sort good chords by their corresponding cycle length
    std::sort(goodChords.begin(), goodChords.end(), [](auto &c1, auto &c2) { return c1.length < c2.length; });

    // loop over chords
    auto coveredNodes = std::vector<bool>(n_, false);
    // add ocw when nodes are not covered yet
    for (auto cChord : goodChords) {
        if (!coveredNodes[cChord.start] && !coveredNodes[cChord.end]) {
            auto newOcw = ocw.splitOnChord(cChord.start, cChord.end, cChord.inner, cChord.crossing);
            assert(newOcw.isSimple());
            result.push_back(newOcw);
            for (uint64_t i = 1; i < newOcw.size() - 1; i++) {
                coveredNodes[newOcw.getIthNode(i)] = true;
            }
        }
    }
    if (result.empty()) {
        return {ocw};
    }
    return result;
}

std::vector<OddClosedWalk> OcwSeparator::getViolatedOcws(NetworKit::node source) {
    currentSource_ = source;
    std::vector<OddClosedWalk> violatedOcws;

    if (blacklist_[currentSource_])
        return violatedOcws;
    runDijsktra();
    if (intermediateVertices_.empty()) {
        return violatedOcws;
    }

    // source not blocklisted and at least one intermediate
    std::vector<OddClosedWalk> startOcws = {};
    for (int i = 0; i < std::ssize(intermediateVertices_) && i < 10; ++i) {
        auto intermediate = intermediateVertices_[i];
        startOcws.emplace_back(intermediateVertexToOcw(intermediate));
    }

    for (const auto &currentStartOcw : startOcws) {
        const std::vector<OddClosedWalk> simpleOcws = simpleOddCycles(currentStartOcw);
        assert(!simpleOcws.empty());
        for (const auto &cSimpleOcw : simpleOcws) {
            const std::vector<OddClosedWalk> extractedOcws = extractCycles(cSimpleOcw);
            assert(!extractedOcws.empty());
            for (auto &extractedOcw : extractedOcws) {
                violatedOcws.push_back(extractedOcw);
            }
        }
    }

    return violatedOcws;
}

void OcwSeparator::blacklistIntReachable(NetworKit::node source) {
    assert(source < originalGraph_.upperNodeIdBound());
    blacklist_[source] = true;

    std::vector<NetworKit::node> dfsStack;
    dfsStack.push_back(source);

    while (!dfsStack.empty()) {
        auto next = dfsStack.back();
        dfsStack.pop_back();

        for (unsigned int i = 0; i < originalGraph_.degree(next); i++) {
            auto [neighbor, edgeId] = originalGraph_.getIthNeighborWithId(next, i);
            if ((weights_[edgeId] == 0.0 || weights_[edgeId] == 1.0) && (!blacklist_[neighbor])) {
                blacklist_[neighbor] = true;
                dfsStack.push_back(neighbor);
            }
        }
    }
}

} // namespace sms
