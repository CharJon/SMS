#include "sms/instance/gridded_ising.hpp"

namespace sms {

gridIsing::gridIsing(const Graph &g, int dim, int gridL, double scalingFactor, double offset)
    : Ising(g, scalingFactor, offset), dim_(dim), gridLength_(gridL) {
    layerSize_ = gridL * gridL; // size of the 2d layer, required for movement in 3d grids
}

torusIsing::torusIsing(const Graph &g, int dim, int gridL, double scalingFactor, double offset)
    : gridIsing(g, dim, gridL, scalingFactor, offset) {
    layerSize_ = gridL * gridL; // size of the 2d layer, required for movement in 3d grids
}

node gridUp(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({l, (r == 0) ? gridLength - 1 : r - 1, c}, gridLength);
}

node gridDown(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({l, (r + 1) % gridLength, c}, gridLength);
}

node gridLeft(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({l, r, (c == 0) ? gridLength - 1 : c - 1}, gridLength);
}

node gridRight(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({l, r, (c + 1) % gridLength}, gridLength);
}

node gridBack(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({(l + 1) % gridLength, r, c}, gridLength);
}

node gridFront(node u, int gridLength) {
    auto [l, r, c] = getCoordinates(u, gridLength);
    return getNodeAt({(l == 0) ? gridLength - 1 : l - 1, r, c}, gridLength);
}

bool gridIsing::isConsistent() {
    auto spins = static_cast<int>(pow(gridLength_, dim_));

    if ((n_ != spins) || (m_ != spins * dim_ - dim_ * pow(gridLength_, dim_ - 1))) {
        return false;
    }

    for (auto u : graph_.nodeRange()) {
        // only check down edge if not in the lowest row
        if (!graph_.hasEdge(u, down(u)) && (getRow(u) != static_cast<unsigned int>(gridLength_ - 1))) {
            return false;
        }

        // only check right edge if not in right column
        if (!graph_.hasEdge(u, right(u)) && (getColumn(u) != static_cast<unsigned int>(gridLength_ - 1))) {
            return false;
        }

        if (dim_ == 3) {
            // only check back edge if not in last layer
            if (!graph_.hasEdge(u, back(u)) && (getLayer(u) != static_cast<unsigned int>(gridLength_ - 1)))
                return false;
        }
    }

    return true;
}

nlohmann::ordered_json gridIsing::getInstanceInformation() {
    nlohmann::ordered_json j;

    j["num spins"] = n_;
    j["num interactions"] = m_;
    j["scaling factor"] = scalingFactor_;
    j["offset"] = offset_;
    j["dimension"] = dim_;
    j["grid length"] = gridLength_;

    return j;
}

void gridIsing::printInstanceInformation(std::ostream &out) {
    out << "---------- Gridded Ising instance statistics -----------------" << std::endl;
    auto j = getInstanceInformation();
    out << j.dump(4) << std::endl;
    out << "--------------------------------------------------------------" << std::endl;
}

bool torusIsing::isConsistent() {
    auto spins = static_cast<int>(pow(gridLength_, dim_));

    if ((n_ != spins) || (m_ != spins * dim_)) {
        return false;
    }

    for (auto u : graph_.nodeRange()) {
        if (!graph_.hasEdge(u, down(u)) || !graph_.hasEdge(u, right(u))) {
            return false;
        }

        if (dim_ == 3) {
            if (!graph_.hasEdge(u, back(u)))
                return false;
        }
    }

    return true;
}

std::optional<std::tuple<int, int>> computeTorusParams(int n, int m) {
    int dim = -1;
    if (m % 2 == 0) {
        std::cout << "Number of Interactions suggests two dimensional torus. Attempting to verify..." << std::endl;
        int gridLength = static_cast<int>(std::sqrt(n));
        dim = 2;
        if ((gridLength * gridLength == n) && (n * dim == m)) {
            std::cout << "Verified torus as two dimensional" << std::endl;
            return {{2, gridLength}};
        } else
            std::cout << "Could not verify torus as two dimensional" << std::endl;
    }

    if (m % 3 == 0) {
        std::cout << "Number of Interactions suggests three dimensional torus. Attempting to verify..." << std::endl;
        int gridLength = static_cast<int>(std::cbrt(n));
        dim = 3;
        if ((gridLength * gridLength * gridLength == n) && (n * dim == m)) {
            std::cout << "Verified torus as three dimensional" << std::endl;
            return {{3, gridLength}};
        }
    }

    std::cout << std::endl;

    if ((dim != -1) && (n * dim != m)) {
        std::cerr << "Number of Interactions does not match computed dimension:" << std::endl;
        std::cerr << "Dimension                                  : " << dim << std::endl;
        std::cerr << "#Spins                                     : " << n << std::endl;
        std::cerr << "#Interactions                              : " << m << std::endl;
        std::cerr << "Expected #Interactions (#Spins * Dimension): " << n * dim << std::endl;
    } else {
        std::cerr << "Could not determine torus dimension for given interactions and spins:" << std::endl;
        std::cerr << "Number of Interactions should be divisible by the dimension of the torus." << std::endl;
        std::cerr << "#Interactions % 2: " << m % 2 << std::endl;
        std::cerr << "#Interactions % 3: " << m % 3 << std::endl;
        std::cerr << "Number of Spins should have an integer root of the dimension of the torus." << std::endl;
        std::cerr << "#Spins ^ 1/2: " << std::sqrt(n) << std::endl;
        std::cerr << "#Spins ^ 1/3: " << std::cbrt(n) << std::endl;
    }

    return {};
}

std::optional<std::tuple<int, int>> computeGridParams(int n, int m) {
    int dim = -1;
    int gridLength = -1;

    gridLength = static_cast<int>(std::sqrt(n));
    if (gridLength * gridLength == n) {
        std::cout << "Number of Spin suggests two dimensional grid. Attempting to verify..." << std::endl;
        dim = 2;
        if (dim * (n - pow(gridLength, dim - 1)) == m) {
            std::cout << "Verified grid as two dimensional." << std::endl;
            return {{dim, gridLength}};
        } else {
            std::cout << "Could not verify grid as two dimensional. Trying three dimensional." << std::endl;
        }
    }

    gridLength = static_cast<int>(std::cbrt(n));
    if (gridLength * gridLength * gridLength == n) {
        std::cout << "Number of Spin suggests three dimensional grid. Attempting to verify..." << std::endl;
        dim = 3;
        if (dim * (n - pow(gridLength, dim - 1)) == m) {
            std::cout << "Verified grid as three dimensional." << std::endl;
            return {{dim, gridLength}};
        } else {
            std::cout << "Could not verify grid as three dimensional." << std::endl;
        }
    }

    if ((dim != -1) && (n * dim != m)) {
        std::cerr << "Number of Interactions does not match computed dimension:" << std::endl;
        std::cerr << "Dimension                                  : " << dim << std::endl;
        std::cerr << "#Spins                                     : " << n << std::endl;
        std::cerr << "#Interactions                              : " << m << std::endl;
        std::cerr << "Expected #Interactions (Dimension * (#Spins - #Spins ^ [(Dimension - 1) / Dimension])): "
                  << (dim * (n - pow(gridLength, dim - 1))) << std::endl;
    } else {
        std::cerr << "Could not determine grid dimension for given interactions and spins:" << std::endl;
        std::cerr << "Number of Interactions should have an integer root of the dimension of the torus." << std::endl;
        std::cerr << "#Spins ^ 1/2: " << std::sqrt(n) << std::endl;
        std::cerr << "#Spins ^ 1/3: " << std::cbrt(n) << std::endl;
    }

    return {};
}

gridCoordinates getCoordinates(node u, int gridLength) {
    return {getLayer(u, gridLength), getRow(u, gridLength), getColumn(u, gridLength)};
}

node getNodeAt(gridCoordinates coords, int gridLength) {
    auto uGridLength = static_cast<node>(gridLength);
    return coords.layer * (uGridLength * uGridLength) + coords.row * uGridLength + coords.column;
}

} // namespace sms
