#ifndef SMS_ODD_CLOSED_WALK_HPP
#define SMS_ODD_CLOSED_WALK_HPP

#include "optional"
#include <vector>

#include "networkit/graph/Graph.hpp"

namespace sms {

struct ParityEdge {
    NetworKit::node u;
    NetworKit::node v;
    bool isStay;
};

/***
 * Saving data for an odd closed walk
 * Start and end node are the same
 ***/
class OddClosedWalk {
private:
    std::vector<NetworKit::node> nodes_;
    std::vector<bool> isStay_;
    uint64_t numCrosses_;

public:
    OddClosedWalk() = delete;

    explicit OddClosedWalk(NetworKit::node startNode, int64_t approxSize = 3) {
        nodes_.reserve(approxSize + 1); // +1 for end node
        isStay_.reserve(approxSize);
        nodes_.push_back(startNode);
        numCrosses_ = 0;
    }

    void addEdge(NetworKit::node nextNode, bool stay) {
        nodes_.push_back(nextNode);
        isStay_.push_back(stay);
        numCrosses_ += 1 - stay;
    }

    // start and end are the same, they are counted only once
    uint64_t size() const { return nodes_.size() - 1; }

    uint64_t numCrosses() const { return numCrosses_; }

    NetworKit::node getIthNode(uint64_t i) const {
        assert(i < nodes_.size());
        return nodes_[i];
    }

    bool ithEdgeIsStay(uint64_t i) const {
        assert(i < isStay_.size());
        return isStay_[i];
    }

    bool isSimple() const {
        std::unordered_set<NetworKit::node> seen;
        if (nodes_.empty()) {
            return false;
        }
        for (unsigned int i = 0; i < nodes_.size() - 1; ++i) {
            auto n = nodes_[i];
            if (seen.contains(n)) {
                return false;
            }
            seen.insert(n);
        }
        return true;
    }

    bool isValid() const {
        bool startAndEndSame = (nodes_[0] == nodes_[size()]);
        bool correctLength = (isStay_.size() == size());
        uint64_t numCross = 0;
        for (auto i : isStay_) {
            numCross += (1 - i);
        }
        return startAndEndSame && correctLength && (numCross == numCrosses_) && (numCrosses_ % 2 == 1);
    }

    bool isValid(const NetworKit::Graph &g) const {
        bool isCycle = true;
        for (uint64_t i = 0; i < size(); i++) {
            isCycle &= g.hasEdge(nodes_[i], nodes_[i + 1]);
        }

        return isCycle & isValid();
    }

    /***
     * Extract an inner cycle from start and end position inclusive
     * @param start
     * @param end
     * @return
     */
    OddClosedWalk extract(uint64_t start, uint64_t end) const {
        assert(nodes_[start] == nodes_[end]);
        OddClosedWalk ocw(nodes_[start]);
        for (uint64_t i = start + 1; i <= end; i++) {
            ocw.addEdge(nodes_[i], isStay_[i - 1]);
        }
        assert(ocw.isValid());
        return ocw;
    }

    /***
     * Splits a cycle by using the given chord
     * @param start
     * @param end
     * @param inner if inner, the result is the cycle start -> start+1 -> ... -> end, else end -> end + 1 -> ... ->
     * start
     * @param cross shall the chord be of crossing type
     * @return
     */
    OddClosedWalk splitOnChord(uint64_t start, uint64_t end, bool inner, bool cross) const;

    void inverse() {
        std::reverse(nodes_.begin(), nodes_.end());
        std::reverse(isStay_.begin(), isStay_.end());
    }
};

std::optional<OddClosedWalk> violatedOddSelectionClosedWalk(const std::vector<NetworKit::node> &nodes,
                                                            const std::vector<double> &weights);

} // namespace sms

#endif // SMS_ODD_CLOSED_WALK_HPP
