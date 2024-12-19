#ifndef SMS_BIPARTITION_HPP
#define SMS_BIPARTITION_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace sms {
/*
 * Partition type
 * Casting bool to Partition is always safe
 * Casting partition to bool is safe, if value is not kUNASSIGNED
 */
enum Partition : uint8_t {
    kZERO = 0,
    kONE = 1,
    kUNASSIGNED = 2,
};

const std::string kPartitionZeroName = "partition_0";
const std::string kPartitionOneName = "partition_1";

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> bipartitionToIndexSets(const std::vector<bool> &bipartition);

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> bipartitionToIndexSets(const std::vector<bool> &bipartition,
                                                                               const std::vector<uint64_t> &indexMap);

} // namespace sms

#endif // SMS_BIPARTITION_HPP
