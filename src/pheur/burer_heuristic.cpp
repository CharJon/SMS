#include <numbers>

#include "networkit/auxiliary/Timer.hpp"

#include "sms/pheur/burer_heurisitic.hpp"

namespace sms {

Burer::Burer(const NetworKit::Graph &g, uint64_t seed)
    : MaxCutSolver(g),
      randomGenerator_(seed),
      uniformMinusPiToPi_(-std::numbers::pi, std::numbers::pi),
      tempAssignments_(graph_.upperNodeIdBound(), false),
      theta_(graph_.upperNodeIdBound(), 0.0),
      assignments_(graph_.upperNodeIdBound(), false),
      runtimeLimit_(std::chrono::duration<double>::max()) {
    assert(g.upperNodeIdBound() == g.numberOfNodes());
}

void Burer::runRepeatedly() {
    t.restart();
    setW1norm();
    bestValue_ = -std::numeric_limits<double>::max();
    for (numOuterLoops_ = 0; keepGoing(); ++numOuterLoops_) {
        initRandomTheta();
        coreRun();
    }
    t.stop();
}

void Burer::run() {
    t.restart();
    setW1norm();
    bestValue_ = -std::numeric_limits<double>::max();
    initRandomTheta();
    coreRun();
    t.stop();
}

void Burer::run(std::vector<double> theta) {
    t.restart();
    setW1norm();
    bestValue_ = -std::numeric_limits<double>::max();
    theta_ = std::move(theta);
    coreRun();
    t.stop();
}

void Burer::run(const std::vector<bool> &start) {
    t.restart();
    setW1norm();
    perturbedToTheta(start);
    coreRun();
    t.stop();
}

void Burer::coreRun() {
    // Algorithm 1 from Section 4
    int numNotImproved = 0;
    while (numNotImproved <= maxNumNonImproving_ && keepGoing()) {
        theta_ = gradientDescent(theta_);
        findBestPartition(theta_);

        if (localSearch_) {
            // TODO: report here before local search, as this might be slow for big graphs
            All1Swap(oneMoveTolerance_);
            All2Swap(twoMoveTolerance_);
        }

        if (improvesOver(bestTempValue_, bestValue_)) {
            newSolution(tempAssignments_, bestTempValue_); // update and report new solution
            numNotImproved = 0;
        } else {
            numNotImproved += 1;
        }

        perturbedToTheta(tempAssignments_);
    }
}

void Burer::initRandomTheta() {
    std::generate(theta_.begin(), theta_.end(),
                  [&]() { return uniformMinusPiToPi_(randomGenerator_) + std::numbers::pi; });
}

void Burer::perturbedToTheta(const std::vector<bool> &assignments) {
    for (const auto u : graph_.nodeRange()) {
        theta_[u] = std::numbers::pi * assignments[u] + perturbationK_ * randomBetweenMinusPiAndPi();
    }
}

void Burer::setW1norm() {
    w1norm_ = 0.0;
    for (const auto e : graph_.edgeWeightRange()) {
        w1norm_ += 2.0 * fabs(e.weight); // Count both directions of edge
    }
}

std::vector<bool> Burer::getBestSolution() const {
    return assignments_;
}

std::vector<double> Burer::gradientDescent(std::vector<double> theta) {
    // *** Setup variables
    // The backtracking alpha value for each iteration's line search
    double bt_alpha = alpha_init_;
    // element-wise cosine of the theta vector
    std::vector<double> cos_theta(graph_.upperNodeIdBound());
    // element-wise sine of the theta vector
    std::vector<double> sin_theta(graph_.upperNodeIdBound());
    // negative of cosine differences incident to each node
    std::vector<double> dH(graph_.upperNodeIdBound());
    // current objective value of nonlinear optimization problem
    double f = loadTheta(theta, cos_theta, sin_theta, dH);

    for (int i = 0; i < max_opt_iter_; ++i) {
        // Compute the current gradient
        std::vector<double> gradient = computeGradient(sin_theta, cos_theta);

        // Stop the gradient descent if the gradient is too close to 0
        double norm_gradient = 0.0;
        for (auto u : graph_.nodeRange()) {
            norm_gradient += gradient[u] * gradient[u];
        }
        if (norm_gradient / w1norm_ < g_tolerance_) {
            break;
        }

        // Determine the descent direction
        std::vector<double> desc = computeDescentDirection(gradient, dH);
        double gradient_times_desc = scalarProduct(gradient, desc);

        // Use backtracking Armijo line-search to determine a good step size
        std::vector<double> new_theta(graph_.upperNodeIdBound());
        int numback;
        double recent_f = -1.0;
        for (numback = 1; numback <= maxback_; ++numback) {
            // Compute the new theta with step size alpha
            for (auto u : graph_.nodeRange()) {
                new_theta[u] = theta[u] + bt_alpha * desc[u];
            }

            // Update the cosine and sine vectors, as well as dH and the objective. Exit on Armijo condition.
            recent_f = loadTheta(new_theta, cos_theta, sin_theta, dH);
            if (recent_f <= f + gamma_ * bt_alpha * gradient_times_desc) {
                break;
            }

            // Exponential scaling on step size
            bt_alpha *= tau_;
        }

        double f_prev = f;
        f = recent_f;
        theta = new_theta; // Copy over the current theta

        // Stop the gradient descent if there was not enough change in the
        // objective value
        double rel_change = fabs(f - f_prev) / (1.0 + fabs(f_prev));
        if (rel_change < f_tolerance_) {
            break;
        }

        // Step size update lifted from circut code; not mentioned in paper...
        if (numback <= 2) {
            if (max_alpha_ < 2.0 * bt_alpha) {
                bt_alpha = max_alpha_;
            } else {
                bt_alpha *= 2.0;
            }
        }
    }
    // Modulo the angles to be between 0 and 2*PI
    modulo2Pi(theta);
    return theta;
}

void Burer::findBestPartition(const std::vector<double> &theta) {
    // add angles to a vector of angle/node pairs. Sort on angle.
    std::vector<std::pair<double, NetworKit::node> > angles;
    for (auto u : graph_.nodeRange()) {
        angles.emplace_back(theta[u], u);
    }
    angles.emplace_back(2 * std::numbers::pi, graph_.upperNodeIdBound());
    std::sort(angles.begin(), angles.end());

    // Determine initial set inclusion, and setup variables to be updated later.
    std::fill(tempAssignments_.begin(), tempAssignments_.end(), false);

    unsigned indexSmallerPi = 0; // iterates angles in [0, pi]
    unsigned indexLargerPi = 0;  // iterates angles with theta > pi

    // set indexLargerPi to first node with theta > pi and make first bi-partition set angles in [0, pi]
    for (unsigned i = 0; i < angles.size(); ++i) {
        if (angles[i].first <= std::numbers::pi) {
            tempAssignments_[angles[i].second] = true; // The first set has angles in [0, pi]
        } else {
            indexLargerPi = i; // indexLargerPi is index of the first node with theta > pi
            break;
        }
    }

    // Fill in diffWeights_ and bestTempValue_ from tempAssignments_, and copy them over to the temporary solution we'll
    // be// updating in the search for a better solution.
    bestTempValue_ = loadAssignment(tempAssignments_);
    double curr_weight = bestTempValue_;
    std::vector<bool> curr_assignments(tempAssignments_);
    std::vector<double> curr_diff_weights(diffWeights_);

    // Compute the optimal cut by exhaustively searching through the possible cuts
    NetworKit::node node;
    if (angles[indexSmallerPi].first <= angles[indexLargerPi].first - std::numbers::pi) {
        node = angles[indexSmallerPi].second; // We will exclude this guy
        ++indexSmallerPi;
    } else {
        node = angles[indexLargerPi].second; // We will include this guy
        ++indexLargerPi;
    }
    while (node != graph_.upperNodeIdBound()) {
        // See if this is a new best cut; if so, copy as the current solution
        updateCutValues(node, curr_assignments, curr_diff_weights, curr_weight);

        if (improvesOver(curr_weight, bestTempValue_)) {
            // update best solution
            tempAssignments_ = curr_assignments;
            // update best value
            bestTempValue_ = curr_weight;
            diffWeights_ = curr_diff_weights;
        }
        // Update the cut angle
        if (angles[indexSmallerPi].first <= angles[indexLargerPi].first - std::numbers::pi) {
            node = angles[indexSmallerPi].second; // We will exclude this guy
            ++indexSmallerPi;
        } else {
            node = angles[indexLargerPi].second; // We will include this guy
            ++indexLargerPi;
        }
    }
}

double Burer::loadAssignment(const std::vector<bool> &assignment) {
    double value = 0.0;
    diffWeights_.assign(graph_.upperNodeIdBound(), 0.0);

    for (const auto e : graph_.edgeWeightRange()) {
        auto u = e.u;
        auto v = e.v;
        auto weight = e.weight;
        if (assignment[u] == assignment[v]) {
            // These two nodes are in the same set
            diffWeights_[u] += weight;
            diffWeights_[v] += weight;
        } else {
            // These two nodes are in different sets
            value += weight;
            diffWeights_[u] -= weight;
            diffWeights_[v] -= weight;
        }
    }
    return value;
}

// Does all the updates required when a node has its cut set inclusion flipped
void Burer::updateCutValues(NetworKit::node node, std::vector<bool> &assignment, std::vector<double> &diff_weights,
                            double &objective) const {
    objective += diff_weights[node];
    assignment[node] = !assignment[node];
    diff_weights[node] = -diff_weights[node];

    // Iterate the set of all neighbors for the node
    for (auto [neighbor, weight] : graph_.weightNeighborRange(node)) {
        if (assignment[node] == assignment[neighbor]) {
            diff_weights[neighbor] += 2.0 * weight;
        } else {
            diff_weights[neighbor] -= 2.0 * weight;
        }
    }
}

void Burer::updateCutValues(NetworKit::node node) {
    updateCutValues(node, tempAssignments_, diffWeights_, bestTempValue_);
}

double Burer::loadTheta(const std::vector<double> &theta, std::vector<double> &cos_theta,
                        std::vector<double> &sin_theta, std::vector<double> &diff) {
    updateAngles(theta, cos_theta, sin_theta);
    return cosineDiff(theta, cos_theta, sin_theta, diff);
}

void Burer::updateAngles(const std::vector<double> &theta, std::vector<double> &cos_theta,
                         std::vector<double> &sin_theta) {
    assert(theta.size() == graph_.upperNodeIdBound());
    assert(cos_theta.size() == graph_.upperNodeIdBound());
    assert(sin_theta.size() == graph_.upperNodeIdBound());
    for (unsigned i = 0; i < theta.size(); i++) {
        cos_theta[i] = fastCos(theta[i]);
        sin_theta[i] = fastSin(theta[i]);
    }
}

double Burer::cosineDiff(const std::vector<double> &theta, const std::vector<double> &cos_theta,
                         const std::vector<double> &sin_theta, std::vector<double> &diff) {
    assert(theta.size() == graph_.upperNodeIdBound());
    assert(cos_theta.size() == graph_.upperNodeIdBound());
    assert(sin_theta.size() == graph_.upperNodeIdBound());

    diff = std::vector<double>(theta.size(), 0.0);
    double objective = 0.0;

    for (const auto e : graph_.edgeWeightRange()) {
        auto u = e.u;
        auto v = e.v;
        auto weight = e.weight;
        double scaled_cos_diff = weight * ((cos_theta[u] * cos_theta[v]) + (sin_theta[u] * sin_theta[v]));
        objective += scaled_cos_diff;
        diff[u] -= scaled_cos_diff;
        diff[v] -= scaled_cos_diff;
    }
    return objective;
}

// Perform all first available 1-swaps.
void Burer::All1Swap(double tolerance) {
    bool move_made = true;
    while (move_made) {
        move_made = false;
        for (auto i : graph_.nodeRange()) {
            if (diffWeights_[i] > tolerance && improvesOverTempSolution(i)) {
                updateCutValues(i);
                move_made = true;
                break;
            }
        }
    }
}

// Perform all first available 2-swaps.
void Burer::All2Swap(double tolerance) {
    bool move_made = true;
    while (move_made) {
        move_made = false;
        for (auto e : graph_.edgeWeightRange()) {
            int u = e.u;
            int v = e.v;
            double weight = e.weight;
            double benefit = diffWeights_[u] + diffWeights_[v]
                             - 2.0 * (tempAssignments_[u] == tempAssignments_[v] ? 1 : -1) * weight;
            if (benefit > tolerance && improvesOverTempSolution(benefit)) {
                updateCutValues(u);
                updateCutValues(v);
                move_made = true;
                break;
            }
        }
    }
}
std::vector<double> Burer::computeDescentDirection(const std::vector<double> &gradient, const std::vector<double> &dH) {
    // Determine the descent direction via scaling; we determined this behavior by actually looking at the circut
    // code as this is not mentioned in the Burer2002 paper
    std::vector<double> desc(graph_.upperNodeIdBound());
    assert(gradient.size() == desc.size() && dH.size() == desc.size());

    double divisor = 1.0;
    for (auto u : graph_.nodeRange()) {
        divisor = std::max(divisor, dH[u]);
    }
    for (auto u : graph_.nodeRange()) {
        desc[u] = -gradient[u] / divisor;
    }
    return desc;
}

std::vector<double> Burer::computeGradient(const std::vector<double> &sin_theta, const std::vector<double> &cos_theta) {
    std::vector<double> gradient(graph_.upperNodeIdBound(), 0.0); // Gradient of f(theta)
    for (auto e : graph_.edgeWeightRange()) {
        auto u = e.u;
        auto v = e.v;
        auto weight = e.weight;
        double scaled_sin_diff = weight * (sin_theta[v] * cos_theta[u] - cos_theta[v] * sin_theta[u]);
        gradient[u] += scaled_sin_diff;
        gradient[v] -= scaled_sin_diff; // sin(theta[u] - theta[v]) = -sin(theta[v] - theta[u])
    }
    return gradient;
}

} // namespace sms