#ifndef SMS_MQLIB_HPP
#define SMS_MQLIB_HPP

#include "mqlib/problem/instance.h"
#include "networkit/graph/Graph.hpp"

std::vector<mqlib::Instance::InstanceTuple> mqlibEdgelist(const NetworKit::Graph &g);

#endif // SMS_MQLIB_HPP
