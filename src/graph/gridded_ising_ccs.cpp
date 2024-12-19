#include "sms/graph/gridded_ising_ccs.hpp"

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/instance/ising.hpp"

namespace sms {

void GriddedIsingChordlessCycles::run() {
    for (auto u : instance.getGraph().nodeRange()) {
        auto fourCCs = getFourCyclesAt(u);
        fourHoles.insert(fourHoles.end(), fourCCs.begin(), fourCCs.end());

        if (instance.getDim() == 3) {
            auto sixCCs = getSixCyclesAt(u);
            sixCycles.insert(sixCycles.end(), sixCCs.begin(), sixCCs.end());
        }

        if (instance.getNumberOfSpins() > 8) {
            auto eightCCs = getEightCyclesAt(u);
            eightCycles.insert(eightCycles.end(), eightCCs.begin(), eightCCs.end());
        }
    }

    run_ = true;
}

bool GriddedIsingChordlessCycles::hasRun() const {
    return run_;
}

// assume the given node is upper left corner the 2d case
// for 3d case only cycle that start with right and back are listed
std::vector<fourcycle> GriddedIsingChordlessCycles::getFourCyclesAt(node u) {
    std::vector<fourcycle> cycles;

    auto u_2 = instance.right(u);
    auto u_3 = instance.down(u_2);
    auto u_4 = instance.left(u_3);

    cycles.push_back({u, u_2, u_3, u_4});

    if (instance.getDim() == 3) {
        u_2 = instance.back(u);
        u_3 = instance.down(u_2);
        u_4 = instance.front(u_3);

        cycles.push_back({u, u_2, u_3, u_4});

        u_2 = instance.back(u);
        u_3 = instance.right(u_2);
        u_4 = instance.front(u_3);

        cycles.push_back({u, u_2, u_3, u_4});
    }

    return cycles;
}

std::vector<sixcycle> GriddedIsingChordlessCycles::getSixCyclesAt(node u) {
    std::vector<sixcycle> cycles;

    auto u_2 = instance.back(u);
    auto u_3 = instance.right(u_2);
    auto u_4 = instance.down(u_3);
    auto u_5 = instance.front(u_4);
    auto u_6 = instance.left(u_5);

    cycles.push_back({u, u_2, u_3, u_4, u_5, u_6});

    u_2 = instance.right(u);
    u_3 = instance.down(u_2);
    u_4 = instance.back(u_3);
    u_5 = instance.left(u_4);
    u_6 = instance.up(u_5);

    cycles.push_back({u, u_2, u_3, u_4, u_5, u_6});

    u_2 = instance.down(u);
    u_3 = instance.back(u_2);
    u_4 = instance.right(u_3);
    u_5 = instance.up(u_4);
    u_6 = instance.front(u_5);

    cycles.push_back({u, u_2, u_3, u_4, u_5, u_6});

    u_2 = instance.back(u);
    u_3 = instance.left(u_2);
    u_4 = instance.down(u_3);
    u_5 = instance.front(u_4);
    u_6 = instance.right(u_5);

    cycles.push_back({u, u_2, u_3, u_4, u_5, u_6});

    return cycles;
}

std::vector<eightcycle> GriddedIsingChordlessCycles::getEightCyclesAt(node u) {
    std::vector<eightcycle> cycles;

    auto u_2 = instance.right(u);
    auto u_3 = instance.right(u_2);
    auto u_4 = instance.down(u_3);
    auto u_5 = instance.down(u_4);
    auto u_6 = instance.left(u_5);
    auto u_7 = instance.left(u_6);
    auto u_8 = instance.up(u_7);

    cycles.push_back({u, u_2, u_3, u_4, u_5, u_6, u_7, u_8});

    if (instance.getDim() == 3) {
        u_2 = instance.back(u);
        u_3 = instance.back(u_2);
        u_4 = instance.down(u_3);
        u_5 = instance.down(u_4);
        u_6 = instance.front(u_5);
        u_7 = instance.front(u_6);
        u_8 = instance.up(u_7);

        cycles.push_back({u, u_2, u_3, u_4, u_5, u_6, u_7, u_8});

        u_2 = instance.back(u);
        u_3 = instance.back(u_2);
        u_4 = instance.right(u_3);
        u_5 = instance.right(u_4);
        u_6 = instance.front(u_5);
        u_7 = instance.front(u_6);
        u_8 = instance.up(u_7);

        cycles.push_back({u, u_2, u_3, u_4, u_5, u_6, u_7, u_8});
    }

    return cycles;
}

} // namespace sms