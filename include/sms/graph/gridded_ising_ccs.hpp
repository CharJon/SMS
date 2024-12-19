#ifndef SMS_GRIDDED_ISING_CCS_HPP
#define SMS_GRIDDED_ISING_CCS_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/instance/gridded_ising.hpp"

namespace sms {
using NetworKit::Graph;
using NetworKit::node;
using fourcycle = std::array<node, 4>;
using sixcycle = std::array<node, 6>;
using eightcycle = std::array<node, 8>;

class GriddedIsingChordlessCycles {
public:
    const torusIsing &instance;

    std::vector<fourcycle> fourHoles;
    std::vector<sixcycle> sixCycles;
    std::vector<eightcycle> eightCycles;

    explicit GriddedIsingChordlessCycles(const torusIsing &instance) : instance(instance) {
        fourHoles = std::vector<fourcycle>();
        sixCycles = std::vector<sixcycle>();
        eightCycles = std::vector<eightcycle>();
    }

    void run();

    bool hasRun() const;

    // get ccs beginning at given node
    std::vector<fourcycle> getFourCyclesAt(node u);

    std::vector<sixcycle> getSixCyclesAt(node u);

    std::vector<eightcycle> getEightCyclesAt(node u);

private:
    bool run_ = false;
};

} // namespace sms

#endif // SMS_GRIDDED_ISING_CCS_HPP
