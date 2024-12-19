#include "sms/instance/convert.hpp"

namespace sms {

MaxCut frustrationIndexToMaxCut(const FrustrationIndex &frustrationIndex) {
    NetworKit::Graph mcGraph(frustrationIndex.getNumberOfNodes(), true);

    double offset = 0.0;
    double scaling = -1.0;

    for (auto e : frustrationIndex.getGraph().edgeWeightRange()) {
        mcGraph.addEdge(e.u, e.v, -1.0 * e.weight);
        if (e.weight == -1.0)
            offset += -1.0;
    }

    auto inst = MaxCut(mcGraph);
    inst.setOffset(offset);
    inst.setScalingFactor(scaling);
    return inst;
}

QUBO maxCutToQUBO(const MaxCut &maxcut) {
    const int dim = maxcut.getNumberOfVertices();
    QUBO q(dim);
    const Graph &g = maxcut.getGraph();

    for (auto e : g.edgeWeightRange()) {
        q.setValue(e.u, e.v, e.weight);
        q.setValue(e.v, e.u, e.weight);
    }

    for (auto u : g.nodeRange()) {
        double sum = 0;
        for (auto neighbor : g.neighborRange(u)) {
            sum -= g.weight(u, neighbor);
        }
        q.setValue(u, u, sum);
    }

    q.setOffset(maxcut.getOffset());
    q.setScalingFactor(-maxcut.getScalingFactor());

    return q;
}

QUBO maxCutToQUBOrooted(const MaxCut &maxcut) {
    const auto dim = maxcut.getNumberOfVertices() - 1;
    QUBO q(dim);
    const Graph &g = maxcut.getGraph();

    for (int u = 0; u < dim; u++) {
        double sum = 0;
        for (auto [neighbor, weight] : g.weightNeighborRange(u)) {
            if (neighbor < static_cast<NetworKit::node>(dim))
                q.setValue(u, neighbor, weight);
            sum -= g.weight(u, neighbor);
        }
        q.setValue(u, u, sum);
    }

    q.setOffset(q.getOffset() + maxcut.getOffset());
    q.setScalingFactor(-maxcut.getScalingFactor());

    return q;
}

MaxCut isingToMaxCut(const Ising &ising) {
    NetworKit::Graph mcGraph(ising.getNumberOfSpins(), true, false);
    double wSum = 0;

    for (auto e : ising.getGraph().edgeWeightRange()) {
        mcGraph.addEdge(e.u, e.v, -2 * e.weight);
        wSum += e.weight;
    }

    auto inst = MaxCut(mcGraph);
    inst.setOffset(wSum);
    return inst;
}

MaxCut QUBOToMaxCut(const QUBO &qubo) {
    NetworKit::node rootNode = qubo.getDim(); // the new node, needed for transform
    NetworKit::Graph g(rootNode + 1, true, false);

    for (int a = 0; a < qubo.getDim(); a++) {
        for (int b = a + 1; b < qubo.getDim(); b++) {
            double w = qubo.getValue(a, b) + qubo.getValue(b, a);
            w = w / 2;
            if (w != 0.0) {
                assert(!g.hasEdge(a, b));
                g.addEdge(a, b, w);
            }
        }
    }

    for (int a = 0; a < qubo.getDim(); a++) {
        auto rootEdgeWeight = -1 * qubo.getValue(a, a) - g.weightedDegree(a);
        if (rootEdgeWeight != 0.0) {
            g.addEdge(a, rootNode, rootEdgeWeight);
        }
    }

    auto inst = MaxCut(g);
    inst.setScalingFactor(-1 * qubo.getScalingFactor());
    inst.setOffset(qubo.getOffset());
    return inst;
}

std::vector<std::array<node, 4>> graphToTorusIsing::getFourHolesAt(node u, node v) {
    std::vector<std::array<node, 4>> res;

    if (graph.hasEdge(u, v)) {
        for (auto uNeighbor : graph.neighborRange(u)) {
            for (auto vNeighbor : graph.neighborRange(v)) {
                if ((uNeighbor != vNeighbor) && (uNeighbor != v) && (vNeighbor != u)
                    && graph.hasEdge(uNeighbor, vNeighbor)) {
                    res.push_back({u, v, vNeighbor, uNeighbor});
                }
            }
        }
    }

    return res;
}

graphToTorusIsing::graphToTorusIsing(const Graph &g) : graph(g) {
    assert(g.numberOfNodes() == g.upperNodeIdBound());

    auto n = g.numberOfNodes();
    auto m = g.numberOfEdges();

    auto params = computeTorusParams(static_cast<int>(n), static_cast<int>(m));

    run_ = false;

    if (!params.has_value()) {
        isGrid_ = false;
    } else {
        auto [d, l] = params.value();
        dim_ = d;
        gridLength_ = l;

        unassigned = 2 + n;
        gridToOrig = std::vector<node>(n, unassigned);
        origToGrid = std::vector<node>(n, unassigned);
    }
}

void graphToTorusIsing::run() {
    run_ = true;

    for (auto v : graph.nodeRange()) {
        if (static_cast<int>(graph.degree(v)) != 2 * dim_) {
            isGrid_ = false;
            return;
        }
    }

    if (dim_ == 2) {

        if (!fixFirstCycle()) // set up the first cycle
        {
            isGrid_ = false;
            return;
        }

        // as the first cycle is now fixed, we can start fixing the first row
        if (!fixFirstRow()) // set up the first row
        {
            isGrid_ = false;
            return;
        }

        // start fixing the rows

        node currentRow = gridDown(0, gridLength_);

        while (currentRow != 0) {
            if (!fixRow(currentRow)) // set up the row
            {
                isGrid_ = false;
                return;
            }

            currentRow = gridDown(currentRow, gridLength_);
        }

        isGrid_ = true;
    }
    if (dim_ == 3) {
        if (!fixFirstCube()) // set up the first cube
        {
            isGrid_ = false;
            return;
        }

        // as the first cube is now fixed, we can start fixing the first row
        if (!fixFirstRow3D()) // set up the first row
        {
            isGrid_ = false;
            return;
        }

        node currentRow = gridDown(0, gridLength_);

        while (currentRow != gridUp(0, gridLength_)) { // end early as last row won't have any unassigned nodes
            if (!fixRow3D(currentRow))                 // set up the row
            {
                isGrid_ = false;
                return;
            }

            currentRow = gridDown(currentRow, gridLength_);
        }

        // check edges from up to down of the last row (these are not yet checked)
        currentRow = 0;
        do {
            if (!graph.hasEdge(gridToOrig[gridUp(currentRow, gridLength_)], gridToOrig[currentRow])
                or !graph.hasEdge(gridToOrig[gridBack(currentRow, gridLength_)],
                                  gridToOrig[gridBack(gridUp(currentRow, gridLength_), gridLength_)])) {
                isGrid_ = false;
                return;
            }
            currentRow = gridRight(currentRow, gridLength_);
        } while (currentRow != 0);

        // row going back

        currentRow = gridBack(0, gridLength_);

        while (currentRow != gridFront(0, gridLength_)) { // end early as last row won't have any unassigned nodes
            if (!fixRowFrontFace3D(currentRow))           // set up the row
            {
                isGrid_ = false;
                return;
            }

            currentRow = gridBack(currentRow, gridLength_);
        }

        // check edges from front to back of the last row (these are not yet checked)
        currentRow = 0;
        do {
            if (!graph.hasEdge(gridToOrig[gridFront(currentRow, gridLength_)], gridToOrig[currentRow])
                or !graph.hasEdge(gridToOrig[gridDown(currentRow, gridLength_)],
                                  gridToOrig[gridDown(gridFront(currentRow, gridLength_), gridLength_)])) {
                isGrid_ = false;
                return;
            }
            currentRow = gridRight(currentRow, gridLength_);
        } while (currentRow != 0);

        // now the remaining nodes
        node currentLayer = gridBack(gridDown(0, gridLength_), gridLength_);

        while (currentLayer != gridBack(gridUp(0, gridLength_), gridLength_)) {
            // one layer
            currentRow = currentLayer;

            while (currentRow
                   != gridFront(currentLayer, gridLength_)) { // end early as last row won't have any unassigned nodes
                if (!fixRowFrontUpperFace3D(currentRow))      // set up the row
                {
                    isGrid_ = false;
                    return;
                }

                currentRow = gridBack(currentRow, gridLength_);
            }

            currentRow = gridDown(0, gridLength_);
            do {
                if (!graph.hasEdge(gridToOrig[gridFront(currentRow, gridLength_)], gridToOrig[currentRow])
                    or !graph.hasEdge(gridToOrig[gridDown(currentRow, gridLength_)],
                                      gridToOrig[gridDown(gridFront(currentRow, gridLength_), gridLength_)])) {
                    isGrid_ = false;
                    return;
                }
                currentRow = gridRight(currentRow, gridLength_);
            } while (currentRow != gridDown(0, gridLength_));

            currentLayer = gridDown(currentLayer, gridLength_);
        }

        // some final edges form the top-most to the bottom-most layer need to be checked still

        node a = 0;
        node b = 0;

        while (a != gridLeft(0, gridLength_)) {
            b = a;
            while (b != gridFront(a, gridLength_)) {
                if (!graph.hasEdge(gridToOrig[gridUp(b, gridLength_)], gridToOrig[b])) {
                    isGrid_ = false;
                    return;
                }
                b = gridBack(b, gridLength_);
            }
            a = gridRight(a, gridLength_);
        }

        isGrid_ = true;
    }
}

torusIsing graphToTorusIsing::getTorusInstance() {
    assert(run_);
    assert(isGrid_);
    if (!isGrid_) {
        throw std::runtime_error("Can not convert to torus instance as graph is not a grid.");
    }

    NetworKit::Graph griddedGraph(graph.numberOfNodes(), true, false);
    for (auto e : graph.edgeWeightRange()) {
        griddedGraph.addEdge(origToGrid[e.u], origToGrid[e.v], e.weight);
    }

    return torusIsing(griddedGraph, dim_, gridLength_);
}

bool graphToTorusIsing::isGrid() const {
    assert(run_);
    return isGrid_;
}

bool graphToTorusIsing::hasRun() const {
    return run_;
}

std::vector<node> graphToTorusIsing::getCommonNeighbours(node u, node v) {
    std::vector<node> res;

    for (auto uN : graph.neighborRange(u)) {
        for (auto vN : graph.neighborRange(v)) {
            if (uN == vN)
                res.push_back(uN);
        }
    }
    return res;
}

// map node 0 and its immediate neighbours
bool graphToTorusIsing::fixFirstCycle() {

    // just set those, does not matter as it is the first node
    gridToOrig[0] = 0;
    origToGrid[0] = 0;

    for (auto u : graph.neighborRange(0)) {
        for (auto v : graph.neighborRange(0)) {
            if (v != u) {
                auto comNs = getCommonNeighbours(u, v);

                if (comNs.size() == 2) // need to find two common neighbours: 0 and one new one
                {
                    auto neigh = comNs[0] == 0 ? comNs[1] : comNs[0]; // get the non-zero common neighbour

                    // as nothing is fixed yet we can set the first cycle arbitrarily
                    gridToOrig[gridDown(0, gridLength_)] = v;  // v will be node below 0
                    gridToOrig[gridRight(0, gridLength_)] = u; // u will be to the right
                    gridToOrig[gridRight(gridDown(0, gridLength_),
                                         gridLength_)] = neigh; // their common neighbour will be down and right from 0

                    std::cout << v << " " << u << " " << neigh << " 0" << std::endl;

                    origToGrid[v] = gridDown(0, gridLength_);
                    origToGrid[u] = gridRight(0, gridLength_);
                    origToGrid[neigh] = gridRight(gridDown(0, gridLength_), gridLength_);
                    return true; // one cycle has been fixed so we are done
                }
            }
        }
    }

    // if we cant even fix the first cycle we are done here, the run function should check this to terminated is
    // possible
    return false;
}

// at this point the first cycle (i.e. 4 nodes should be fixed), starting at the cycle we fix the first row
bool graphToTorusIsing::fixFirstRow() {

    node upperLeft = 1; // upper left corner of the current cycle to fix (in grid space)

    while (upperLeft != 0) {
        // lower left is also fixed
        auto lowerLeft = gridDown(upperLeft, gridLength_);

        // look at the existing ccs for the edges from upper left to its down neighbour
        // both should already be fixed
        auto realCCs = getFourHolesAt(gridToOrig[upperLeft], gridToOrig[lowerLeft]);

        // for 2d each edge is part of 2, for 3d of 4 4-holes, return false if otherwise
        if (realCCs.size() != static_cast<unsigned int>(2 * (dim_ - 1))) {
            return false;
        }

        // one of the two real ccs will be already completely known, the other will contain two unassigned vertices
        std::array<node, 4> cycle{};
        // currently inefficient, but should work
        for (auto cc : realCCs) {
            for (auto v : cc) {
                if (origToGrid[v] == unassigned)
                    cycle = cc;
            }
        }

        // after we found the unassigned cycle we can assign it
        for (auto v : cycle) {
            if (origToGrid[v] == unassigned) {
                if (graph.hasEdge(
                        v,
                        gridToOrig[upperLeft])) // if the vertex is connected to upperLeft it is right of upperLeft
                {
                    gridToOrig[gridRight(upperLeft, gridLength_)] = v;
                    origToGrid[v] = gridRight(upperLeft, gridLength_);
                }

                if (graph.hasEdge(
                        v,
                        gridToOrig[lowerLeft])) // if the vertex is connected to lowerLeft it is right of lowerLeft
                {
                    gridToOrig[gridRight(lowerLeft, gridLength_)] = v;
                    origToGrid[v] = gridRight(lowerLeft, gridLength_);
                }
            }
        }

        // set upperLeft to the next node
        upperLeft = gridRight(upperLeft, gridLength_);
    }

    // if we end here all is good
    return true;
}

// fix a row, assuming that upperLeft and its right neighbour are already fixed (i.e. the first row is fixed)
bool graphToTorusIsing::fixRow(node initUpperLeft) {
    node upperLeft = initUpperLeft; // upper left corner of the current cycle to fix (in grid space)

    do {
        // upper right is also fixed
        auto upperRight = gridRight(upperLeft, gridLength_);

        // look at the existing ccs for the edges from upper left to its down neighbour
        // both should already be fixed
        auto realCCs = getFourHolesAt(gridToOrig[upperLeft], gridToOrig[upperRight]);

        // for 2d each edge is part of 2, for 3d of 4 4-holes, return false if otherwise
        if (realCCs.size() != static_cast<unsigned int>(2 * (dim_ - 1))) {
            return false;
        }

        // one of the two real ccs will be already completely known, the other will contain one or two unassigned
        // vertices
        std::array<node, 4> cycle{};
        // currently inefficient, but should work
        for (auto cc : realCCs) {
            for (auto v : cc) {
                if (origToGrid[v] == unassigned)
                    cycle = cc;
            }
        }

        // after we found the unassigned cycle we can assign it
        for (auto v : cycle) {
            if (origToGrid[v] == unassigned) {
                if (graph.hasEdge(
                        v,
                        gridToOrig[upperLeft])) // if the vertex is connected to upperLeft it is right of upperLeft
                {
                    gridToOrig[gridDown(upperLeft, gridLength_)] = v;
                    origToGrid[v] = gridDown(upperLeft, gridLength_);
                }

                if (graph.hasEdge(
                        v,
                        gridToOrig[upperRight])) // if the vertex is connected to lowerLeft it is right of lowerLeft
                {
                    gridToOrig[gridDown(upperRight, gridLength_)] = v;
                    origToGrid[v] = gridDown(upperRight, gridLength_);
                }
            }
        }

        // set upperLeft to upperRight
        upperLeft = upperRight;
    } while (upperLeft != initUpperLeft);

    // if we end here all is good
    return true;
}

// fix the first 2 by 2 by 2 cube for the 3d grid case
bool graphToTorusIsing::fixFirstCube() {
    // just set those, does not matter as it is the first node
    gridToOrig[0] = 0;
    origToGrid[0] = 0;

    bool done = false;

    for (auto u : graph.neighborRange(0)) {
        for (auto v : graph.neighborRange(0)) {
            if ((v != u) && !done) {
                auto comNs = getCommonNeighbours(u, v);

                if (comNs.size() == 2) // need to find two common neighbours: 0 and one new one
                {
                    auto neigh = comNs[0] == 0 ? comNs[1] : comNs[0]; // get the non-zero common neighbour

                    // as nothing is fixed yet we can set the first cycle arbitrarily
                    gridToOrig[gridDown(0, gridLength_)] = v;  // v will be node below 0
                    gridToOrig[gridRight(0, gridLength_)] = u; // u will be to the right
                    gridToOrig[gridRight(gridDown(0, gridLength_),
                                         gridLength_)] = neigh; // their common neighbour will be down and right from 0

                    origToGrid[v] = gridDown(0, gridLength_);
                    origToGrid[u] = gridRight(0, gridLength_);
                    origToGrid[neigh] = gridRight(gridDown(0, gridLength_), gridLength_);
                    // one cycle has been fixed so we are done
                    done = true;
                }
            }
        }
    }

    /* now that we fixed the first cycle, we can try to expand it to the first cube
     * for better understanding lets name all those nodes (labelled u,v,w,x clockwise starting from 0):
     * u -- v           0 -- 1
     * |    |    <=>    |    |
     * x -- w           l -- l+1        (l = gridLength)
     */

    node u = 0;
    node v = gridToOrig[gridRight(0, gridLength_)];
    node w = gridToOrig[gridDown(gridRight(0, gridLength_), gridLength_)];
    node x = gridToOrig[gridDown(0, gridLength_)];

    /* we assume that this cycle now if the front face of cube, now lets find the upper face
     * to do so we look at the 4-holes which include (u, v).
     * one will contain {u,v,w,x}, the three other will only include {u,v}
     * from those three one will be the upper face of the cube above our current cube
     * the two remaining can both be labelled as the upper face, we can identify by checking that the two new node in
     * those cycles have each one common nighbour with x and w, those nighbours are then the last nodes of the cube
     * (back(x) and (back(w))
     *
     * In general we will call the 4 nodes of the back face y,z,a,b:
     * y -- z           back(u) -- back(v)
     * |    |    <=>      |          |
     * b -- a           back(x) -- back(w)
     *
     */

    // step one get all candidate cycles
    auto candidates = getFourHolesAt(u, v);

    if (candidates.size() != 4) // wrong number of cycles means no grid
    {
        return false;
    }

    for (auto c : candidates) {
        std::vector<node> newNodes;
        for (auto cNode : c) {
            if (origToGrid[cNode] == unassigned) // means that node is not assigned
            {
                newNodes.push_back(
                    cNode); // new_nodes include the two nodes for the cycle which have not been assigned yet
            }
        }

        if (newNodes.size() == 2) // possible candidate for the upper face
        {
            node y = graph.hasEdge(newNodes[0], u) ? newNodes[0] : newNodes[1];
            node z = graph.hasEdge(newNodes[0], v) ? newNodes[0] : newNodes[1];
            node a = unassigned;
            node b = unassigned;

            auto comXY = getCommonNeighbours(x, y); // this will be node b, u
            auto comWZ = getCommonNeighbours(w, z); // this will be node a, v

            if (comWZ.size() == 2 && comXY.size() == 2) // check that the sizes work out
            {
                // now fix a and b
                b = comXY[0] == u ? comXY[1] : comXY[0];
                a = comWZ[0] == v ? comWZ[1] : comWZ[0];

                // finally fix the back face
                assert(comWZ.size() == 2);
                assert(comXY.size() == 2);
                gridToOrig[gridBack(origToGrid[u], gridLength_)] = y;
                gridToOrig[gridBack(origToGrid[v], gridLength_)] = z;
                gridToOrig[gridBack(origToGrid[w], gridLength_)] = a;
                gridToOrig[gridBack(origToGrid[x], gridLength_)] = b;

                origToGrid[y] = gridBack(origToGrid[u], gridLength_);
                origToGrid[z] = gridBack(origToGrid[v], gridLength_);
                origToGrid[a] = gridBack(origToGrid[w], gridLength_);
                origToGrid[b] = gridBack(origToGrid[x], gridLength_);

                return true;
            }
        }
    }

    // if we cant even fix the first cube we are done here, the run function should check this to terminated is possible
    return false;
}

bool graphToTorusIsing::fixFirstRow3D() {
    node upperLeft = 1; // upper left corner of the current cube to fix (in grid space)

    while (upperLeft != gridLeft(0, gridLength_)) {
        // lower left is also fixed
        auto lowerLeft = gridDown(upperLeft, gridLength_);
        // as well as the back nodes for lower and upper left
        auto backUpperLeft = gridBack(upperLeft, gridLength_);
        auto backLowerLeft = gridBack(lowerLeft, gridLength_);

        /* w.l.o.g we assume that the left face is the one already fixed, so we now need to the right face of this cube
         * this can be done by the ccs of the upper egde (upperLeft, backUpperLeft) and then checking for common
         * neighbours of the two new nodes of this 4-hole (similar to fixing the first cube)
         */

        // look at the existing ccs for the edges from upper left to its back neighbour
        // both should already be fixed, there should be 4 cycles, 2 which are totally fixed
        auto candidates = getFourHolesAt(gridToOrig[upperLeft], gridToOrig[backUpperLeft]);

        // two of the cycles will be completely fixed, two can have 2 un fixed nodes
        if (candidates.size() != 4) // wrong number of cycles means no grid
        {
            return false;
        }
        bool success = false;
        for (auto c : candidates) {
            if (!success) {
                std::vector<node> new_nodes;
                for (auto c_node : c) {
                    if (origToGrid[c_node] == unassigned) // means that node is not assigned
                    {
                        new_nodes.push_back(
                            c_node); // new_nodes include the two nodes for the cycle which have not been assigned yet
                    }
                }
                if (new_nodes.size()
                    == 2) // possible candidate for the upper face, should have the two upper right nodes as new nodes
                {
                    node upperRight = graph.hasEdge(new_nodes[0], gridToOrig[upperLeft]) ? new_nodes[0] : new_nodes[1];
                    node backUpperRight =
                        graph.hasEdge(new_nodes[0], gridToOrig[backUpperLeft]) ? new_nodes[0] : new_nodes[1];
                    node lowerRight = unassigned;
                    node backLowerRight = unassigned;

                    auto com_upperRight_lowerLeft =
                        getCommonNeighbours(upperRight,
                                            gridToOrig[lowerLeft]); // this will be lowerRight
                    auto com_backUpperRight_backLowerLeft =
                        getCommonNeighbours(backUpperRight,
                                            gridToOrig[backLowerLeft]); // this will be node backLowerRight

                    if (com_upperRight_lowerLeft.size() == 2
                        && com_backUpperRight_backLowerLeft.size() == 2) // check that the sizes work out
                    {
                        // now fix a and b
                        lowerRight = com_upperRight_lowerLeft[0] == gridToOrig[upperLeft] ? com_upperRight_lowerLeft[1]
                                                                                          : com_upperRight_lowerLeft[0];
                        backLowerRight = com_backUpperRight_backLowerLeft[0] == gridToOrig[backUpperLeft]
                                             ? com_backUpperRight_backLowerLeft[1]
                                             : com_backUpperRight_backLowerLeft[0];

                        // finally fix the back face

                        gridToOrig[gridRight(upperLeft, gridLength_)] = upperRight;
                        gridToOrig[gridRight(backUpperLeft, gridLength_)] = backUpperRight;
                        gridToOrig[gridRight(lowerLeft, gridLength_)] = lowerRight;
                        gridToOrig[gridRight(backLowerLeft, gridLength_)] = backLowerRight;

                        origToGrid[upperRight] = gridRight(upperLeft, gridLength_);
                        origToGrid[backUpperRight] = gridRight(backUpperLeft, gridLength_);
                        origToGrid[lowerRight] = gridRight(lowerLeft, gridLength_);
                        origToGrid[backLowerRight] = gridRight(backLowerLeft, gridLength_);

                        success = true;
                    }
                }
            }
        }
        if (success) {
            upperLeft = gridRight(upperLeft, gridLength_);
        } else {
            return false;
        }
    }

    // now all should be assigned, check that the last edges exist
    auto lowerLeft = gridDown(upperLeft, gridLength_);
    // as well as the back nodes for lower and upper left
    auto backUpperLeft = gridBack(upperLeft, gridLength_);
    auto backLowerLeft = gridBack(lowerLeft, gridLength_);

    bool closing_edges = graph.hasEdge(gridToOrig[upperLeft], gridToOrig[gridRight(upperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[lowerLeft], gridToOrig[gridRight(lowerLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backUpperLeft], gridToOrig[gridRight(backUpperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backLowerLeft], gridToOrig[gridRight(backLowerLeft, gridLength_)]);

    // if we end here all is good
    return closing_edges;
}

bool graphToTorusIsing::fixRow3D(node initUpperLeft) {
    node upperLeft = initUpperLeft; // upper left corner of the current cube to fix (in grid space)

    /* Should be similar to fixing the first row. Except that we suppose the upper face of the cube is fully fixed,
     * so we need to construct the lower one.
     * */

    while (upperLeft != gridLeft(initUpperLeft, gridLength_)) {
        // back upper left is also fixed
        auto backUpperLeft = gridBack(upperLeft, gridLength_);
        // as well as the right nodes of the uer face
        auto backUpperRight = gridRight(backUpperLeft, gridLength_);
        auto upperRight = gridRight(upperLeft, gridLength_);

        // look at the existing ccs for the edges from upper left to its back neighbour
        // both should already be fixed, there should be 4 cycles, 3 which are totally fixed (as the first row is fixed)
        auto candidates = getFourHolesAt(gridToOrig[upperLeft], gridToOrig[backUpperLeft]);

        // two of the cycles will be completely fixed, one should have 2 unfixed nodes
        if (candidates.size() != 4) // wrong number of cycles means no grid
        {
            return false;
        }
        bool success = false;
        for (auto c : candidates) {
            if (!success) {
                std::vector<node> new_nodes;
                for (auto c_node : c) {
                    if (origToGrid[c_node] == unassigned) // means that node is not assigned
                    {
                        new_nodes.push_back(
                            c_node); // new_nodes include the two nodes for the cycle which have not been assigned yet
                    }
                }
                if ((new_nodes.size() == 2)
                    or (upperLeft != initUpperLeft)) // only possible candidate for the left face, should have the two
                                                     // lower left nodes as new nodes, however for every cube after the
                                                     // first one no new nodes will be found this way
                {
                    node lowerLeft = unassigned;
                    node backLowerLeft = unassigned;
                    if (upperLeft == initUpperLeft) { // required at first cube of the row
                        lowerLeft = graph.hasEdge(new_nodes[0], gridToOrig[upperLeft]) ? new_nodes[0] : new_nodes[1];
                        backLowerLeft =
                            graph.hasEdge(new_nodes[0], gridToOrig[backUpperLeft]) ? new_nodes[0] : new_nodes[1];
                    } else // required later on
                    {
                        lowerLeft = gridToOrig[gridDown(upperLeft, gridLength_)];
                        backLowerLeft = gridToOrig[gridDown(backUpperLeft, gridLength_)];
                    }

                    node lowerRight = unassigned;
                    node backLowerRight = unassigned;

                    auto com_upperRight_lowerLeft = getCommonNeighbours(gridToOrig[upperRight],
                                                                        lowerLeft); // this will be lowerRight
                    auto com_backUpperRight_backLowerLeft =
                        getCommonNeighbours(gridToOrig[backUpperRight],
                                            backLowerLeft); // this will be node backLowerRight

                    if (com_upperRight_lowerLeft.size() == 2
                        && com_backUpperRight_backLowerLeft.size() == 2) // check that the sizes work out
                    {
                        // now fix a and b
                        lowerRight = com_upperRight_lowerLeft[0] == gridToOrig[upperLeft] ? com_upperRight_lowerLeft[1]
                                                                                          : com_upperRight_lowerLeft[0];
                        backLowerRight = com_backUpperRight_backLowerLeft[0] == gridToOrig[backUpperLeft]
                                             ? com_backUpperRight_backLowerLeft[1]
                                             : com_backUpperRight_backLowerLeft[0];

                        // finally fix the lower face
                        if (upperLeft == initUpperLeft) // also fix lower lefts node
                        {
                            gridToOrig[gridDown(upperLeft, gridLength_)] = lowerLeft;
                            gridToOrig[gridDown(backUpperLeft, gridLength_)] = backLowerLeft;

                            origToGrid[lowerLeft] = gridDown(upperLeft, gridLength_);
                            origToGrid[backLowerLeft] = gridDown(backUpperLeft, gridLength_);
                        }

                        gridToOrig[gridDown(upperRight, gridLength_)] = lowerRight;
                        gridToOrig[gridDown(backUpperRight, gridLength_)] = backLowerRight;

                        origToGrid[lowerRight] = gridDown(upperRight, gridLength_);
                        origToGrid[backLowerRight] = gridDown(backUpperRight, gridLength_);

                        success = true;
                    } else {
                        return false; // as we only get on candidate we can return false, if it does not work out
                    }
                }
            }
        }
        if (success) {
            upperLeft = gridRight(upperLeft, gridLength_); // move further in the row
        } else {
            return false;
        }
    }

    // now all should be assigned, check that the last edges exist
    auto lowerLeft = gridDown(upperLeft, gridLength_);
    // as well as the back nodes for loewer and upper left
    auto backUpperLeft = gridBack(upperLeft, gridLength_);
    auto backLowerLeft = gridBack(lowerLeft, gridLength_);

    bool closing_edges = graph.hasEdge(gridToOrig[upperLeft], gridToOrig[gridRight(upperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[lowerLeft], gridToOrig[gridRight(lowerLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backUpperLeft], gridToOrig[gridRight(backUpperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backLowerLeft], gridToOrig[gridRight(backLowerLeft, gridLength_)]);

    // if we end here all is good
    return closing_edges;
}

// similar to fix row, but now assumes the front face to be already given indtaed of the upper face
bool graphToTorusIsing::fixRowFrontFace3D(node initUpperLeft) {
    node upperLeft = initUpperLeft; // upper left corner of the current cube to fix (in grid space)

    /* Should be similar to fixing the rows. Except that we suppose the front face of the cube is fully fixed,
     * so we need to construct the back one.
     * */

    while (upperLeft != gridLeft(initUpperLeft, gridLength_)) {
        // right upper, and lower left and right are also fixed
        auto upperRight = gridRight(upperLeft, gridLength_);
        // as well as the right nodes of the uer face
        auto lowerRight = gridDown(upperRight, gridLength_);
        auto lowerLeft = gridDown(upperLeft, gridLength_);

        // look at the existing ccs for the edges from upper left to its back neighbour
        // both should already be fixed, there should be 4 cycles, 3 which are totally fixed (as the first row is fixed)
        auto candidates = getFourHolesAt(gridToOrig[upperLeft], gridToOrig[lowerLeft]);

        // two of the cycles will be completely fixed, one should have 2 unfixed nodes
        if (candidates.size() != 4) // wrong number of cycles means no grid
        {
            return false;
        }
        bool success = false;
        for (auto c : candidates) {
            if (!success) {
                std::vector<node> new_nodes;
                for (auto c_node : c) {
                    if (origToGrid[c_node] == unassigned) // means that node is not assigned
                    {
                        new_nodes.push_back(
                            c_node); // new_nodes include the two nodes for the cycle which have not been assigned yet
                    }
                }
                if ((new_nodes.size() == 2)
                    or (upperLeft != initUpperLeft)) // only possible candidate for the right face, should have the two
                                                     // lower left nodes as new nodes, hoever for every cube after the
                                                     // first one no new nodes will be found this way
                {
                    node backUpperLeft = unassigned;
                    node backLowerLeft = unassigned;
                    if (upperLeft == initUpperLeft) { // required at first cube of the row
                        backUpperLeft =
                            graph.hasEdge(new_nodes[0], gridToOrig[upperLeft]) ? new_nodes[0] : new_nodes[1];
                        backLowerLeft =
                            graph.hasEdge(new_nodes[0], gridToOrig[lowerLeft]) ? new_nodes[0] : new_nodes[1];
                    } else // required later on
                    {
                        backUpperLeft = gridToOrig[gridBack(upperLeft, gridLength_)];
                        backLowerLeft = gridToOrig[gridBack(lowerLeft, gridLength_)];
                    }

                    node backUpperRight = unassigned;
                    node backLowerRight = unassigned;

                    auto com_upperRight_backUpperLeft =
                        getCommonNeighbours(gridToOrig[upperRight],
                                            backUpperLeft); // this will be backUpperRight
                    auto com_lowerRight_backLowerLeft =
                        getCommonNeighbours(gridToOrig[lowerRight],
                                            backLowerLeft); // this will be node backLowerRight

                    if (com_upperRight_backUpperLeft.size() == 2
                        && com_lowerRight_backLowerLeft.size() == 2) // check that the sizes work out
                    {
                        // now fix a and b
                        backUpperRight = com_upperRight_backUpperLeft[0] == gridToOrig[upperLeft]
                                             ? com_upperRight_backUpperLeft[1]
                                             : com_upperRight_backUpperLeft[0];
                        backLowerRight = com_lowerRight_backLowerLeft[0] == gridToOrig[lowerLeft]
                                             ? com_lowerRight_backLowerLeft[1]
                                             : com_lowerRight_backLowerLeft[0];

                        // finally fix the lower face
                        if (upperLeft == initUpperLeft) // back right nodes
                        {
                            gridToOrig[gridBack(upperLeft, gridLength_)] = backUpperLeft;
                            gridToOrig[gridBack(lowerLeft, gridLength_)] = backLowerLeft;

                            origToGrid[backUpperLeft] = gridBack(upperLeft, gridLength_);
                            origToGrid[backLowerLeft] = gridBack(lowerLeft, gridLength_);
                        }

                        gridToOrig[gridBack(upperRight, gridLength_)] = backUpperRight;
                        gridToOrig[gridBack(lowerRight, gridLength_)] = backLowerRight;

                        origToGrid[backUpperRight] = gridBack(upperRight, gridLength_);
                        origToGrid[backLowerRight] = gridBack(lowerRight, gridLength_);

                        success = true;
                    } else {
                        return false; // as we only get on candidate we can return false, if it does not work out
                    }
                }
            }
        }
        if (success) {
            upperLeft = gridRight(upperLeft, gridLength_); // move further in the row
        } else {
            return false;
        }
    }

    // now all should be assigned, check that the last edges exist
    auto lowerLeft = gridDown(upperLeft, gridLength_);
    // as well as the back nodes for loewer and upper left
    auto backUpperLeft = gridBack(upperLeft, gridLength_);
    auto backLowerLeft = gridBack(lowerLeft, gridLength_);

    bool closing_edges = graph.hasEdge(gridToOrig[upperLeft], gridToOrig[gridRight(upperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[lowerLeft], gridToOrig[gridRight(lowerLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backUpperLeft], gridToOrig[gridRight(backUpperLeft, gridLength_)]);

    closing_edges =
        closing_edges && graph.hasEdge(gridToOrig[backLowerLeft], gridToOrig[gridRight(backLowerLeft, gridLength_)]);

    // if we end here all is good
    return closing_edges;
}

// same stuff, now front and upper face are expected to be fully fixed
bool graphToTorusIsing::fixRowFrontUpperFace3D(node initUpperLeft) {
    node upperLeft = initUpperLeft; // upper left corner of the current cube to fix (in grid space)

    /* Should be similar to fixing the rows. Except that we suppose the front face and upper face of the cube is fully
     * fixed, so we need to construct the lower back nodes only
     * */

    while (upperLeft != gridLeft(initUpperLeft, gridLength_)) {
        // right upper, and lower left and right are also fixed
        auto upperRight = gridRight(upperLeft, gridLength_);
        // as well as the right nodes of the uer face
        auto lowerRight = gridDown(upperRight, gridLength_);
        auto lowerLeft = gridDown(upperLeft, gridLength_);

        auto backUpperLeft = gridBack(upperLeft, gridLength_);
        auto backUpperRight = gridBack(upperRight, gridLength_);

        node backLowerLeft = unassigned; // only missing nodes
        node backLowerRight = unassigned;

        auto comLowerLeftBackUpperLeft = getCommonNeighbours(gridToOrig[lowerLeft],
                                                             gridToOrig[backUpperLeft]); // this will be backLowerLeft
        auto comLowerRightBackUpperRight =
            getCommonNeighbours(gridToOrig[lowerRight],
                                gridToOrig[backUpperRight]); // this will be node backLowerRight

        if (comLowerLeftBackUpperLeft.size() == 2
            && comLowerRightBackUpperRight.size() == 2) // check that the sizes work out
        {
            // now fix a and b
            backLowerLeft = comLowerLeftBackUpperLeft[0] == gridToOrig[upperLeft] ? comLowerLeftBackUpperLeft[1]
                                                                                  : comLowerLeftBackUpperLeft[0];
            backLowerRight = comLowerRightBackUpperRight[0] == gridToOrig[upperRight] ? comLowerRightBackUpperRight[1]
                                                                                      : comLowerRightBackUpperRight[0];

            // finally fix the lower face
            gridToOrig[gridBack(lowerLeft, gridLength_)] = backLowerLeft;
            gridToOrig[gridBack(lowerRight, gridLength_)] = backLowerRight;

            origToGrid[backLowerLeft] = gridBack(lowerLeft, gridLength_);
            origToGrid[backLowerRight] = gridBack(lowerRight, gridLength_);

            upperLeft = gridRight(upperLeft, gridLength_); // move further in the row
        } else {
            return false; // as we only get on candidate we can return false, if it does not work out
        }
    }

    // now all should be assigned, check that the last edges exist
    auto lowerLeft = gridDown(upperLeft, gridLength_);
    // as well as the back nodes for loewer and upper left
    auto backUpperLeft = gridBack(upperLeft, gridLength_);
    auto backLowerLeft = gridBack(lowerLeft, gridLength_);

    bool closingEdges = graph.hasEdge(gridToOrig[upperLeft], gridToOrig[gridRight(upperLeft, gridLength_)]);

    closingEdges = closingEdges && graph.hasEdge(gridToOrig[lowerLeft], gridToOrig[gridRight(lowerLeft, gridLength_)]);

    closingEdges =
        closingEdges && graph.hasEdge(gridToOrig[backUpperLeft], gridToOrig[gridRight(backUpperLeft, gridLength_)]);

    closingEdges =
        closingEdges && graph.hasEdge(gridToOrig[backLowerLeft], gridToOrig[gridRight(backLowerLeft, gridLength_)]);

    // if we end here all is good
    return closingEdges;
}

} // namespace sms