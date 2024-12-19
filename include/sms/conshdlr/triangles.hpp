#ifndef SMS_TRIANGLES_HPP
#define SMS_TRIANGLES_HPP

#include <optional>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/odd_closed_walk.hpp"

namespace sms {

class TriangleSeparator {
private:
    NetworKit::Graph lpWeightedGraph_;

    std::vector<std::array<NetworKit::node, 3>> triangles_;

public:
    explicit TriangleSeparator(NetworKit::Graph const *const originalGraph)
        : lpWeightedGraph_(originalGraph->numberOfNodes(), true, false, false) {}

    void updateWeights(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        // assert(u < originalGraph_->numberOfNodes());
        // assert(v < originalGraph_->numberOfNodes());

        // Guarantee numerical stability at 0 and 1
        NetworKit::edgeweight fixedVal = std::min(std::max(0., w), 1.);

        lpWeightedGraph_.setWeight(u, v, fixedVal);
    }

    void addTriangle(NetworKit::node u, NetworKit::node v, NetworKit::node w) { triangles_.push_back({u, v, w}); }

    std::vector<OddClosedWalk> getViolatedTriangles() {
        auto violatedTriangles = std::vector<OddClosedWalk>();

        for (auto t : triangles_) {
            auto [u, v, w] = t;

            std::vector<NetworKit::edgeweight> values = {
                lpWeightedGraph_.weight(u, v),
                lpWeightedGraph_.weight(v, w),
                lpWeightedGraph_.weight(w, u),
            };

            auto violatedWalk = violatedOddSelectionClosedWalk({u, v, w}, values);
            if (violatedWalk.has_value()) {
                violatedTriangles.push_back(violatedWalk.value());
            }
        }

        return violatedTriangles;
    }
};

} // namespace sms

#endif // SMS_TRIANGLES_HPP
