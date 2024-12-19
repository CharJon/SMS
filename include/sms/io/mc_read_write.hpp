#ifndef SMS_MC_READ_WRITE_HPP
#define SMS_MC_READ_WRITE_HPP

#include "networkit/graph/Graph.hpp"

namespace sms {

void writeToMcFile(const NetworKit::Graph &g, const std::string &path);

} // namespace sms

#endif // SMS_MC_READ_WRITE_HPP
