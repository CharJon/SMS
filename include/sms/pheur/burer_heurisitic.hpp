#ifndef SMS_BURER_HEURISITIC_HPP
#define SMS_BURER_HEURISITIC_HPP

#include <chrono>

#include "sms/auxiliary/chrono_io.hpp"
#include "sms/auxiliary/math.hpp"
#include "sms/solver/abstract_solver.hpp"

namespace sms {

class Burer : MaxCutSolver {
public:
    explicit Burer(const NetworKit::Graph &g, uint64_t seed);

    // Any weighted undirected graph is a valid input
    bool applicable() override { return true; }

    void run() override;

    void run(std::vector<double>);

    void run(const std::vector<bool> &start);

    void runRepeatedly();
    void initRandomTheta();

    std::vector<bool> getBestSolution() const override;

    // Always false: Solver is heuristic
    bool optimalityProven() const override { return false; }

    void setSeed(const uint64_t seed) { randomGenerator_.seed(seed); }

    void setTimelimit(const std::chrono::duration<double> time) { runtimeLimit_ = time; }

    void setSolutionValueLimit(const double value) { solutionValueLimit_ = value; }

    const std::vector<double> &getSolutionValues() const { return solutionValues_; }

    const std::vector<double> &getSolutionTimeStamps() const { return solutionTimeStamps_; }

private:
    // randomness
    std::mt19937_64 randomGenerator_;                           // 64 bit random number generator
    std::uniform_real_distribution<double> uniformMinusPiToPi_; // -pi inclusive pi exclusive

    double bestTempValue_ = 0.;
    std::vector<bool> tempAssignments_;
    std::vector<double> theta_;
    std::vector<bool> assignments_;
    std::vector<double> diffWeights_; // Amount you would gain from switching sets / flipping variable
    double w1norm_ = 0.;

    Timer t;
    std::chrono::duration<double> runtimeLimit_;

    std::vector<double> solutionValues_;
    std::vector<double> solutionTimeStamps_;
    double solutionValueLimit_{std::numeric_limits<double>::max()};

    // config
    // Perturbation in the range [-perturbation*PI, perturbation*PI]
    static constexpr double perturbationK_ = 0.2;
    // Stopping value for relative change in f(theta) when optimizing f
    static constexpr double f_tolerance_ = 1e-4;
    // Stopping value for normalized value of g(theta) when optimizing f
    static constexpr double g_tolerance_ = 1e-4;
    // Maximum number of backtracks before termination
    static constexpr int maxback_ = 10;
    // Backtrack multiplier
    static constexpr double tau_ = 0.5;
    // The "beta" value in the Armijo line-search (required fraction of 1st order reduction)
    // see http://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf slide 9
    static constexpr double gamma_ = 0.01;
    // The initial step size in the line search
    static constexpr double alpha_init_ = 1.0;
    // The maximum step size in the line search
    static constexpr double max_alpha_ = 4.0;
    // Maximum number of gradient descent steps when optimizing f(theta)
    static constexpr int max_opt_iter_ = 200;

    // Number of permitted non-improving perturbations to optimal theta before search is stopped.
    static constexpr int maxNumNonImproving_ = 10;
    // Whether we perform greedy 1- and 2-moves after
    static constexpr bool localSearch_ = true;
    // Amount of improvement required to do a 1-move
    static constexpr double oneMoveTolerance_ = 0.01;
    // Amount of improvement required to do a 2-move
    static constexpr double twoMoveTolerance_ = 0.1;

    // stats
    int numOuterLoops_ = 0;

    void coreRun();

    bool keepGoing() const {
        // time limit
        if (t.elapsed() > runtimeLimit_) {
            return false;
        }
        // solution value
        if (bestValue_ >= solutionValueLimit_) {
            return false;
        }
        return true;
    };

    void newSolution(const std::vector<bool> &solution, const double value) {
        if (value > bestValue_) { // if new best -> add to history
            solutionTimeStamps_.push_back(t.elapsed().count());
            solutionValues_.push_back(value);
        }
        // update best solution
        assignments_ = solution;
        // update best value
        bestValue_ = value;
    };

    // Minimize f(theta) using gradient descent with backtracking Armijo line-search. w1norm is the 1-norm of
    // the matrix of weights.
    std::vector<double> gradientDescent(std::vector<double> theta);

    // Find the best bi-partition based on theta
    void findBestPartition(const std::vector<double> &theta);

    // Compute descent direction (desc) and return
    std::vector<double> computeDescentDirection(const std::vector<double> &gradient, const std::vector<double> &dH);

    std::vector<double> computeGradient(const std::vector<double> &, const std::vector<double> &);

    // Given a (new) theta vector, compute the element-wise cosine and sine
    // Return weighted sum of the cosines of the angle differences between paired nodes
    double loadTheta(const std::vector<double> &theta, std::vector<double> &cos_theta, std::vector<double> &sin_theta,
                     std::vector<double> &diff);

    // compute the element-wise cosine and sine
    void updateAngles(const std::vector<double> &theta, std::vector<double> &cos_theta, std::vector<double> &sin_theta);

    // Given a (new) theta vector, compute the weighted sum of the cosines of the angle differences between paired nodes
    double cosineDiff(const std::vector<double> &theta, const std::vector<double> &cos_theta,
                      const std::vector<double> &sin_theta, std::vector<double> &diff);

    // This function initializes diff_weights_ from the given assignment, returns cut value
    double loadAssignment(const std::vector<bool> &assignment);

    // Switch the set of node, updating, the assignments (x) and the diff_weights. Returns the change in cut value.
    void updateCutValues(NetworKit::node node, std::vector<bool> &assignment, std::vector<double> &diff_weights,
                         double &objective) const;

    void updateCutValues(NetworKit::node node);

    bool improvesOver(double weight1, double weight2) const {
        // MightDo : Replace this with function from our math.hpp
        return weight1 - weight2 > 1e-6 && (weight2 <= 0.0 || (weight1 / weight2 - 1.0) >= 1e-10);
    }

    // Identify if a move sufficiently improves the quality of this solution to be taken.
    bool improvesOverTempSolution(NetworKit::node node) const {
        return improvesOver(bestTempValue_ + diffWeights_[node], bestTempValue_);
    }

    bool improvesOverTempSolution(double diff_weight) const {
        return improvesOver(bestTempValue_ + diff_weight, bestTempValue_);
    }

    double randomBetweenMinusPiAndPi() { return uniformMinusPiToPi_(randomGenerator_); }

    void All1Swap(double tolerance);
    void All2Swap(double tolerance);

    // Take bi-partition, map it to angular [-pi, pi] representation and perturb
    // Result is stored in theta_
    void perturbedToTheta(const std::vector<bool> &assignments);

    // Compute 1-norm of vec(W), the vectorization of the weight matrix.
    void setW1norm();

    // Modulo the angles to be between 0 and 2*PI
    void modulo2Pi(std::vector<double> &theta) {
        for (auto u : graph_.nodeRange()) {
            theta[u] -= 2 * std::numbers::pi * floor(theta[u] / (2 * std::numbers::pi)); // Modulo to [0, 2*pi)
        }
    }

    double scalarProduct(const std::vector<double> &v, const std::vector<double> &w) {
        assert(v.size() == w.size());
        double scalarProduct = 0.;
        for (unsigned i = 0; i < v.size(); ++i) {
            scalarProduct += v[i] * w[i];
        }
        return scalarProduct;
    }
};
} // namespace sms
#endif // SMS_BURER_HEURISITIC_HPP
