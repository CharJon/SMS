#ifndef SMS_CONVERT_HPP
#define SMS_CONVERT_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/graph/graphs.hpp"
#include "sms/graph/gridded_ising_ccs.hpp"
#include "sms/instance/frustration_index.hpp"
#include "sms/instance/ising.hpp"
#include "sms/instance/maxcut.hpp"
#include "sms/instance/qubo.hpp"

/***
 * Functions to convert between different instance types.
 * All instances have scaling and offset.
 * A negative scaling factor indicates the problem originally had a different obj-sense than the current one.
 * E.g. if a MaxCut instance has scaling -1, the original problem was a minimization problem.
 * If a QUBO has scaling -1, the original problem was a maximization problem.
 */

namespace sms {

/***
 * Convert a MaxCut instance with n vertices to a QUBO instance with dimension n.
 * Note:
 *   - MaxCut is maximization problem, QUBO is minimization problem, the scaling will be multiplied by -1.
 *   - This will result in a QUBO with at least one optimality-symmetry.
 * @param maxcut
 */
QUBO maxCutToQUBO(const MaxCut &maxcut);

/***
 * Convert a MaxCut instance with n vertices to a QUBO instance with dimension n-1.
 * Note:
 *   - MaxCut is maximization problem, QUBO is minimization problem, the scaling will be multiplied by -1.
 *   - This transformation treats the node with highest id as root node, to avoid adding symmetry to the QUBO.
 * @param maxcut
 */
QUBO maxCutToQUBOrooted(const MaxCut &maxcut);

/***
 * Convert a QUBO instance with dimension d to a MaxCut instance with d+1 vertices.
 * The vertex with the highest id is the newly added node.
 * @param qubo
 */
MaxCut QUBOToMaxCut(const QUBO &qubo);

MaxCut isingToMaxCut(const Ising &ising);

MaxCut frustrationIndexToMaxCut(const FrustrationIndex &frustrationIndex);

class graphToTorusIsing {
public:
    Graph const &graph;
    std::vector<node> gridToOrig; // input grid node, output node in original
    std::vector<node> origToGrid; // input original node, output node in grid

    explicit graphToTorusIsing(const Graph &g);

    void run();

    torusIsing getTorusInstance();

    bool isGrid() const;

    bool hasRun() const;

private:
    bool isGrid_;
    bool run_;
    int dim_;
    int gridLength_;
    node unassigned;

    bool fixFirstCycle();

    bool fixFirstRow();

    bool fixRow(node initUpperLeft);

    // for 3d?
    bool fixFirstCube();

    bool fixFirstRow3D();

    bool fixRow3D(node initUpperLeft);

    bool fixRowFrontFace3D(node initUpperLeft);

    bool fixRowFrontUpperFace3D(node initUpperLeft);

    std::vector<std::array<node, 4>> getFourHolesAt(node u, node v);

    std::vector<node> getCommonNeighbours(node u, node v);
};

} // namespace sms

#endif // SMS_CONVERT_HPP
