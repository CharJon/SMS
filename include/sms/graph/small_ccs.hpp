#ifndef SMS_SMALL_CCS_HPP
#define SMS_SMALL_CCS_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

namespace sms {

using triangle = std::array<NetworKit::node, 3>;
using fourcycle = std::array<NetworKit::node, 4>;

class SmallChordlessCycles {
public:
    NetworKit::Graph const &graph;
    std::vector<triangle> triangles;
    std::vector<fourcycle> fourHoles;

    explicit SmallChordlessCycles(const NetworKit::Graph &graph) : graph(graph) {
        triangles = std::vector<triangle>();
        fourHoles = std::vector<fourcycle>();
    }

    void run();

    bool hasRun() const;

private:
    bool run_ = false;

    void computeSmallCCs();
};

// sort a given triangle according to Lemma 2 of https://arxiv.org/pdf/1309.1051.pdf
triangle sortTriangle(const std::vector<NetworKit::node> &);

// sort a given 4 cycle according to Lemma 2 of https://arxiv.org/pdf/1309.1051.pdf
fourcycle sortFourCycle(const std::vector<NetworKit::node> &, const NetworKit::Graph &);

} // namespace sms

#endif // SMS_SMALL_CCS_HPP
