#ifndef SMS_GRAPHS_STATISTICS_HPP
#define SMS_GRAPHS_STATISTICS_HPP

#include <optional>

#include "nlohmann/json.hpp"

#include "networkit/clique/MaximalCliques.hpp"
#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

#include "sms/auxiliary/statistics.hpp"
#include "sms/graph/small_ccs.hpp"

namespace sms {

nlohmann::ordered_json graphStatsJson(const NetworKit::Graph &g, bool fast = true);

void outputGraphStats(const NetworKit::Graph &g, std::ostream &out, bool fast);

/*
 * Returns the degeneracy of the graph, or std::nullopt if the degeneracy is larger than stopIfLarger
 * Setting stopIfLarger to maxDegree(g) will always return the degeneracy
 */
std::optional<int> degeneracy(const NetworKit::Graph &g, unsigned int stopIfLarger);

} // namespace sms

#endif // SMS_GRAPHS_STATISTICS_HPP
