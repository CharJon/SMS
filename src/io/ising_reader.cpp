#include "sms/io/ising_reader.hpp"

namespace sms {

SgParser::SgParser(const std::string &path, double scaling, double offset) {
    scaling_ = scaling;
    offset_ = offset;

    graph_ = mcFileToGraph(path);
}

Ising SgParser::getInstance() {
    return Ising(graph_, scaling_, offset_);
}

GsgParser::GsgParser(const std::string &path, double scaling, double offset) {
    scaling_ = scaling;
    offset_ = offset;
    graph_ = mcFileToGraph(path);
    torus_ = false;

    isValid();
}

std::optional<torusIsing> GsgParser::getTorus() const {
    if (torus_) {
        return torusIsing(graph_, dim_, gridLength_, scaling_, offset_);
    } else {
        return std::nullopt;
    }
}

std::optional<gridIsing> GsgParser::getGrid() const {
    if (!torus_) {
        return gridIsing(graph_, dim_, gridLength_, scaling_, offset_);
    } else {
        return std::nullopt;
    }
}

void GsgParser::isValid() {
    int n = static_cast<int>(graph_.numberOfNodes());
    int m = static_cast<int>(graph_.numberOfEdges());

    auto torus_params = computeTorusParams(n, m);
    auto grid_params = computeGridParams(n, m);

    if (!torus_params.has_value() && !grid_params.has_value()) {
        throw std::runtime_error("Grid sized could not be determined for given number of spins and interactions");
    }

    if (torus_params.has_value()) {

        auto [d, l] = torus_params.value();
        dim_ = d;
        gridLength_ = l;

        if (n * dim_ != m) {
            throw std::runtime_error("Number of interactions does not match expected value: Expected "
                                     + std::to_string(n * dim_) + " got " + std::to_string(m));
        }

        for (auto i : graph_.nodeRange()) {
            for (int dim = 0; dim < dim_; dim++) {
                node g_tod = static_cast<node>(pow(gridLength_, dim));
                node g_todplus_one = static_cast<node>(pow(gridLength_, dim + 1));
                node u = i - 1;
                node requested = u + g_tod;
                unsigned int depth = u / g_todplus_one + 1;
                node v = requested - ((requested) >= g_todplus_one * depth ? g_todplus_one : 0);

                if (!graph_.hasEdge(u + 1, v + 1)) {
                    throw std::runtime_error("Interaction between " + std::to_string(u + 1) + " and "
                                             + std::to_string(v + 1) + " not specified");
                }
            }
        }

        torus_ = true;
    } else if (grid_params.has_value()) {

        auto [d, l] = grid_params.value();
        dim_ = d;
        gridLength_ = l;

        if ((n - pow(gridLength_, dim_ - 1)) * dim_ != m) {
            throw std::runtime_error("Number of interactions does not match expected value: Expected "
                                     + std::to_string((n - pow(gridLength_, dim_ - 1)) * dim_) + " got "
                                     + std::to_string(m));
        }

        for (auto i : graph_.nodeRange()) {
            for (int dim = 0; dim < dim_; dim++) {
                node g_tod = static_cast<node>(pow(gridLength_, dim));
                node g_todplus_one = static_cast<node>(pow(gridLength_, dim + 1));
                node u = i - 1;
                node requested = u + g_tod;
                unsigned int depth = u / g_todplus_one + 1;
                node v = requested - ((requested) >= g_todplus_one * depth ? g_todplus_one : 0);

                // ignore edges that would close the boundary
                if (u < v) {
                    if (!graph_.hasEdge(u + 1, v + 1)) {
                        throw std::runtime_error("Interaction between " + std::to_string(u + 1) + " and "
                                                 + std::to_string(v + 1) + " not specified");
                    }
                }
            }
        }

        torus_ = false;
    }
}

bool GsgParser::isTorus() const {
    return torus_;
}

} // namespace sms
