#ifndef SMS_HOLES_HPP
#define SMS_HOLES_HPP

#include <optional>

#include "networkit/graph/Graph.hpp"

#include "sms/graph//odd_closed_walk.hpp"

namespace sms {

class HoleSeparator {
private:
    NetworKit::Graph lpWeightedGraph_;

    std::vector<std::array<NetworKit::node, 4>> fourHoles_;

public:
    explicit HoleSeparator(NetworKit::Graph const *const originalGraph)
        : lpWeightedGraph_(originalGraph->numberOfNodes(), true, false, false) {}

    void updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        // assert(u < originalGraph_->numberOfNodes());
        // assert(v < originalGraph_->numberOfNodes());

        // Try to guarantee numerical stability at 0 and 1
        NetworKit::edgeweight fixedVal = std::min(std::max(0., w), 1.);

        lpWeightedGraph_.setWeight(u, v, fixedVal);
    }

    void addFourHole(NetworKit::node a, NetworKit::node b, NetworKit::node c, NetworKit::node d) {
        fourHoles_.push_back({a, b, c, d});
    }

    int numberOfHoles() { return std::ssize(fourHoles_); }

    std::vector<OddClosedWalk> getViolatedFourHoles() {
        auto violatedFourHoles = std::vector<OddClosedWalk>();

        for (auto h : fourHoles_) {
            auto [a, b, c, d] = h;

            std::vector<NetworKit::edgeweight> values = {lpWeightedGraph_.weight(a, b), lpWeightedGraph_.weight(b, c),
                                                         lpWeightedGraph_.weight(c, d), lpWeightedGraph_.weight(a, d)};

            auto violatedWalk = violatedOddSelectionClosedWalk({a, b, c, d}, values);
            if (violatedWalk.has_value()) {
                violatedFourHoles.push_back(violatedWalk.value());
            }
        }

        return violatedFourHoles;
    }
};
} // namespace sms

#endif // SMS_HOLES_HPP
