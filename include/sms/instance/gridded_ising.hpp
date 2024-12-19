#ifndef SMS_GRIDDED_ISING_HPP
#define SMS_GRIDDED_ISING_HPP

#include <optional>
#include <tuple>
#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"
#include "sms/instance/ising.hpp"

namespace sms {
using NetworKit::Edge;
using NetworKit::Graph;
using NetworKit::node;

/* the six functions below return the index of the node located in the indicated direction relative to the given node.
 * The functions go along the edges of the grid. Assume the graphs are 0 based (i.e. node range is 0 to n - 1)*/
node gridUp(node, int);

node gridDown(node, int);

node gridLeft(node, int);

node gridRight(node, int);

node gridBack(node, int); // only 3d

node gridFront(node, int); // only 3d

// return the layer in which node u lies (range from 0 to gridLength - 1) always 0 if dim = 2
inline node getLayer(node u, int gridLength) {
    return u / (gridLength * gridLength);
}

// return row in which u lies within its layer (range from 0 to gridLength - 1)
inline node getRow(node u, int gridLength) {
    return (u - getLayer(u, gridLength) * (gridLength * gridLength)) / gridLength;
}

// return column in which u lies within its layer (range from 0 to gridLength - 1)
inline node getColumn(node u, int gridLength) {
    return u % gridLength;
}

gridCoordinates getCoordinates(node u, int gridLength);

node getNodeAt(gridCoordinates coords, int gridLength);

class gridIsing : public Ising {

public:
    explicit gridIsing(const Graph &g, int dim, int gridL, double scalingFactor = 1, double offset = 0);

    int getDim() const { return dim_; }

    int getGridLength() const { return gridLength_; }

    node up(node u) const { return sms::gridUp(u, gridLength_); }

    node down(node u) const { return sms::gridDown(u, gridLength_); }

    node left(node u) const { return sms::gridLeft(u, gridLength_); }

    node right(node u) const { return sms::gridRight(u, gridLength_); }

    node back(node u) const { // only 3d
        assert(dim_ == 3);
        return sms::gridBack(u, gridLength_);
    }

    node front(node u) const { // only 3d
        assert(dim_ == 3);
        return sms::gridFront(u, gridLength_);
    }

    gridCoordinates getCoordinates(node u) const { return sms::getCoordinates(u, gridLength_); }

    node getNodeAt(gridCoordinates coords) const { return sms::getNodeAt(coords, gridLength_); }

    // check if given graph is consistent with dimension and gridlength
    virtual bool isConsistent();

    nlohmann::ordered_json getInstanceInformation();

    void printInstanceInformation(std::ostream &out);

protected:
    int dim_;
    int gridLength_;
    int layerSize_;

    // return the layer in which node u lies (range from 0 to gridLength - 1) always 0 if dim = 2
    inline node getLayer(node u) const { return sms::getLayer(u, gridLength_); }

    // return row in which u lies within its layer (range from 0 to gridLength - 1)
    inline node getRow(node u) const { return sms::getRow(u, gridLength_); }

    // return column in which u lies within its layer (range from 0 to gridLength - 1)
    inline node getColumn(node u) const { return sms::getColumn(u, gridLength_); }
};

class torusIsing : public gridIsing {
public:
    explicit torusIsing(const Graph &g, int dim, int gridL, double scalingFactor = 1, double offset = 0);

    // check if given graph is consistent with dimension and gridlength
    bool isConsistent() override;
};

/* Given number of nodes n and number of edges m, return a tuple (dim, gridLength) corresponding to the torus
 * described by the given parameters, if possible.
 * */
std::optional<std::tuple<int, int>> computeTorusParams(int n, int m);

/* Given number of nodes n and number of edges m, return a tuple (dim, gridLength) corresponding to the torus
 * described by the given parameters, if possible.
 * */
std::optional<std::tuple<int, int>> computeGridParams(int n, int m);

} // namespace sms
#endif // SMS_GRIDDED_ISING_HPP
