#include "sms/io/mc_read_write.hpp"

namespace sms {

void writeToMcFile(const NetworKit::Graph &g, const std::string &path) {
    std::ofstream file(path);
    if (!file.is_open())
        throw std::runtime_error("File could not be opened: " + path + ".");

    if (g.numberOfSelfLoops() != 0)
        throw std::runtime_error("Graph has selfloops.");
    if (!g.checkConsistency())
        throw std::runtime_error("Graph is inconsistent.");

    file << g.upperNodeIdBound() << " " << g.numberOfEdges() << std::endl;
    for (auto e : g.edgeWeightRange()) {
        assert(e.u != e.v);
        if (e.weight != 0) {
            // Networkit uses 0-based indices, mc uses 1 based
            file << (e.v + 1) << " " << (e.u + 1) << " " << e.weight << std::endl;
        }
    }

    file.close();
}

} // namespace sms