#ifndef SMS_MCSOLUTION_HPP
#define SMS_MCSOLUTION_HPP

#include "nlohmann/json.hpp"

#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/bipartition.hpp"

namespace sms {

NetworKit::edgeweight solutionValue(const NetworKit::Graph &g, const std::vector<bool> &part);

/*
 * Care: Does assert all nodes are partitioned!
 */
NetworKit::edgeweight solutionValue(const NetworKit::Graph &g, const std::vector<Partition> &part);

std::vector<NetworKit::node> allIds(const std::vector<bool> &values, bool match);

nlohmann::json partitionToJson(const NetworKit::Graph &g, const std::vector<bool> &part);

nlohmann::json partitionToJson(const std::vector<NetworKit::node> &existingNodes, const std::vector<bool> &part);

std::vector<bool> partitionFromJson(const nlohmann::json &j);

} // namespace sms

#endif // SMS_MCSOLUTION_HPP
