#include "sms/auxiliary/bipartition.hpp"

namespace sms {

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> bipartitionToIndexSets(const std::vector<bool> &bipartition) {
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> indexSets;
    for (uint64_t i = 0; i < bipartition.size(); i++) {
        auto &p = bipartition[i] ? indexSets.second : indexSets.first;
        p.push_back(i);
    }
    return indexSets;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> bipartitionToIndexSets(const std::vector<bool> &bipartition,
                                                                               const std::vector<uint64_t> &indexMap) {
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> indexSets;
    for (uint64_t i = 0; i < bipartition.size(); i++) {
        auto &p = bipartition[i] ? indexSets.second : indexSets.first;
        p.push_back(indexMap[i]);
    }
    return indexSets;
}

} // namespace sms