#include <cstdlib>
#include <iostream>

#include "networkit/auxiliary/Parallelism.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/io/EdgeListReader.hpp"
#include "networkit/io/MatrixMarketReader.hpp"

#include "sms/auxiliary/cl_parser.hpp"
#include "sms/graph/graphs_statistics.hpp"
#include "sms/instance/convert.hpp"
#include "sms/instance/frustration_index.hpp"
#include "sms/instance/maxcut.hpp"
#include "sms/instance/mc_solution.hpp"
#include "sms/instance/qubo.hpp"
#include "sms/io/io.hpp"
#include "sms/io/ising_reader.hpp"
#include "sms/io/json.hpp"
#include "sms/io/qplib_parser.hpp"
#include "sms/pheur/burer_heurisitic.hpp"
#include "sms/presol/data_reducer_mc.hpp"
#include "sms/presol/presolver_mc.hpp"
#include "sms/solver/basic_ocw.hpp"
#include "sms/solver/heuristic_solver.hpp"

#ifdef SMS_GUROBI
#    include "sms/solver/gurobi_quadratic_solver.hpp"
#endif

void solveInstanceCore(const sms::MaxCut &maxCutInstance, const CLParser &clArgs, nlohmann::ordered_json &stats) {
    auto remainingTime = std::chrono::duration_cast<std::chrono::milliseconds>(clArgs.getTimeLimit());

    bool optimalSolutionFound = true;
    auto presolver = sms::PresolverMC(maxCutInstance.getGraph(), clArgs.getPresolveEmphasis());
    presolver.setTriangleEmphasis(clArgs.getTriangleEmphasis());
    if (clArgs.getOldPresolve())
        presolver.setOldPresolve();
    presolver.run();
    std::cout << "Presolving done, took " << presolver.elapsed() << "." << std::endl;
    std::cout << "Presolver decomposed the graph into " << presolver.numberOfSubgraphs() << " subgraphs." << std::endl;
    if (clArgs.getWarmStartSolutionFile()) {
        std::ifstream warmStartFile(clArgs.getWarmStartSolutionFile().value());
        auto warmStartSolutionJson = nlohmann::json::parse(warmStartFile);
        auto warmStartSolution = sms::partitionFromJson(warmStartSolutionJson);
        // translate solution from original graph to compact graph of instance
        auto transformedSolution = std::vector<bool>(maxCutInstance.getNumberOfVertices(), true);
        for (int i = 0; i < std::ssize(transformedSolution); i++) {
            transformedSolution[i] = warmStartSolution[maxCutInstance.getOriginalNode(i)];
        }
        std::cout << "Reading in warm start solution with value: "
                  << sms::solutionValue(maxCutInstance.getGraph(), transformedSolution) << std::endl;
        presolver.addWarmstartSolution(transformedSolution);
    }

    nlohmann::ordered_json presolverStats = presolver.getCompactStats();
    presolverStats["total runtime [s]"] = std::chrono::duration<double>(presolver.elapsed()).count();
    presolverStats["#subgraphs"] = presolver.numberOfSubgraphs();
    presolverStats["#subgraphs non-trivial"] = presolver.numberOfRemainingSubgraphs();
    presolverStats["#remaining vertices"] = presolver.numRemainingVertices();
    presolverStats["#remaining edges"] = presolver.numRemainingEdges();
    stats["presolver stats"] = presolverStats;
    std::cout << presolverStats.dump(4) << std::endl;

    remainingTime -= presolver.elapsed();
    if (remainingTime <= std::chrono::milliseconds(0)) {
        remainingTime = std::chrono::milliseconds(0);
    }

    nlohmann::ordered_json solverStats;
    Aux::Timer solverTimer;
    // Solving stage
    solverTimer.start();
    for (unsigned int i = 0; i < presolver.numberOfSubgraphs(); i++) {
        // solve each one individually
        const auto &currentReducedSubgraph = presolver.reducedSubgraph(i);
        auto partialSolution = std::vector<sms::Partition>(currentReducedSubgraph.upperNodeIdBound(), sms::kUNASSIGNED);
        auto currentInstance = sms::MaxCut(currentReducedSubgraph);
        std::vector<bool> bestSolution;

        if (currentReducedSubgraph.numberOfNodes() > 1) {
            std::cout << "Solving graph number " << i << " with " << presolver.reducedSubgraph(i).numberOfNodes()
                      << " vertices and " << presolver.reducedSubgraph(i).numberOfEdges() << " edges." << std::endl;

            // TODO: check if special case solver works

            // solve current instance
            if (clArgs.getSolver() == "sms") {
                sms::BasicOcwSolver solver(currentInstance.getGraph());
                solver.setParamFile(clArgs.getParamFile());
                solver.setNodelimit(clArgs.getNodeLimit());
                solver.setIntegral(maxCutInstance.allSolutionsHaveIntegerWeight());
                solver.setSeed(clArgs.getSeed());
                solver.setTimelimit(std::chrono::duration<double>(remainingTime).count());
                if (presolver.hasWarmstartSolution())
                    solver.addWarmStartSolution(presolver.getSubgraphWarmstartSolution(i));
                solver.run();
                optimalSolutionFound &= solver.optimalityProven();
                solverStats["solver " + std::to_string(i)] = solver.getStats();
                bestSolution = solver.getBestSolution();
                // TODO (JC): Bug, remaining time needs to be updated
            } else if (clArgs.getSolver() == "gurobi") {
#ifdef SMS_GUROBI
                std::cout << "Solving with gurobi " << std::endl;
                // solve each one individually
                sms::GurobiQuadraticSolver solver(currentInstance.getGraph());
                solver.setTimelimit(std::chrono::duration<double>(remainingTime).count());
                solver.run();
                solverStats["solver " + std::to_string(i)] = solver.getStats();
                optimalSolutionFound &= solver.optimalityProven();
                bestSolution = solver.getBestSolution();
#else
                std::cerr << "Project not compiled with gurobi flag!" << std::endl;
                exit(EXIT_FAILURE);
#endif
            } else if (clArgs.getSolver() == "heuristic") {
                auto totalNumEdges = static_cast<double>(presolver.numRemainingEdges());
                auto subGraphNumEdges = static_cast<double>(currentInstance.getGraph().numberOfEdges());
                double fraction = subGraphNumEdges / totalNumEdges;
                auto totalRemainingTime = std::chrono::duration_cast<std::chrono::duration<double>>(remainingTime);
                auto runtimeLimit = fraction * totalRemainingTime;
                auto solver = sms::Burer(currentInstance.getGraph(), clArgs.getSeed());
                solver.setTimelimit(runtimeLimit);
                solver.runRepeatedly();

                bestSolution = solver.getBestSolution();
                optimalSolutionFound &= solver.optimalityProven();
                solverStats["solver " + std::to_string(i)]["solution values"] = solver.getSolutionValues();
                solverStats["solver " + std::to_string(i)]["solution timestamps [s]"] = solver.getSolutionTimeStamps();
            } else if (clArgs.getSolver() == "none") {
                optimalSolutionFound = false;
                bestSolution = std::vector<bool>(currentInstance.getNumberOfVertices(), false);
            } else {
                std::cerr << "Unknown solver: " << clArgs.getSolver() << "." << std::endl;
                exit(EXIT_FAILURE);
            }

            // post process results
            for (auto u : currentInstance.getGraph().nodeRange()) {
                bool side = bestSolution[u];
                partialSolution[currentInstance.getOriginalNode(u)] = static_cast<sms::Partition>(side);
            }
        }
        presolver.addSolution(i, partialSolution);
    }
    solverTimer.stop();

    NetworKit::edgeweight bestSolutionValue;
    auto fullSolution = presolver.recoverSolution();
    if (clArgs.getSolutionFile()) {
        // translate solution from compact graph of instance to original graph
        auto transformedSolution = std::vector<bool>(maxCutInstance.maxOriginalNodeId() + 1, false);
        for (int i = 0; i < std::ssize(transformedSolution); i++) {
            if (maxCutInstance.isOriginalNode(i))
                transformedSolution[i] = fullSolution[maxCutInstance.getNewNode(i)];
        }
        auto solJson = sms::partitionToJson(maxCutInstance.originalNodes(), transformedSolution);
        std::ofstream solFile(clArgs.getSolutionFile().value());
        solFile << solJson.dump(4) << std::endl;
    }
    bestSolutionValue = sms::solutionValue(maxCutInstance.getGraph(), fullSolution);
    bestSolutionValue *= maxCutInstance.getScalingFactor();

    solverStats["optimality proven"] = optimalSolutionFound;
    solverStats["runtime [s]"] =
        std::chrono::duration_cast<std::chrono::duration<double>>(solverTimer.elapsed()).count();
    solverStats["best solution value"] = bestSolutionValue;

    stats["solver stats"] = solverStats;
}

void printCommandLineFlags(int argc, char **argv, std::ostream &out = std::cout) {
    out << "Num args: " << argc << ". Args: ";
    for (int i = 0; i < argc; i++) {
        out << argv[i] << " ";
    }
    out << std::endl;
}

void printStartMessage() {
    std::cout << "********************************************" << std::endl;
    std::cout << "*    SMS: A SCIP based MaxCut Solver       *" << std::endl;
    std::cout << "********************************************" << std::endl << std::endl;
}

sms::MaxCut readInstance(const CLParser &clParser) {
    if (!clParser.getFileName().has_value()) {
        std::cout << "No graph as input specified. Will use test graph." << std::endl;
        std::cout << "Call help for more information." << std::endl;
        auto graph = test_graph();
        return sms::MaxCut(graph);
    }

    const auto fileFormat = clParser.getFileFormat();

    // Start with file formats specific to a problem type
    if ((fileFormat == "qplib") && (clParser.getInstType() == "bq")) {
        // ToDo: Seed not used!
        auto qubo = sms::QPLibParser(clParser.getFileName().value()).getInstance();
        return sms::QUBOToMaxCut(qubo);
    }

    if ((fileFormat == "bq") && (clParser.getInstType() == "bq")) {
        // ToDo: Seed not used!
        auto qubo = sms::bqFileToQubo(clParser.getFileName().value());
        return sms::QUBOToMaxCut(qubo);
    }

    if ((fileFormat == "sg") && (clParser.getInstType() == "ising")) {
        // ToDo: Seed not used!
        sms::SgParser parser(clParser.getFileName().value());
        return sms::isingToMaxCut(parser.getInstance());
    }

    if (fileFormat == "gsg") {
        // ToDo: Seed not used!
        sms::GsgParser parser(clParser.getFileName().value());
        // ToDo: Distinguish between torus and grid
        std::cerr << "Not implemented yet." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Input is in a general graph format
    NetworKit::Graph graph;
    if (fileFormat == "mc") {
        graph = mcFileToGraph(clParser.getFileName().value());
    } else if (fileFormat == "wel") {
        NetworKit::EdgeListReader el(' ', 0);
        graph = el.read(clParser.getFileName().value());
    } else if (fileFormat == "wel1") {
        NetworKit::EdgeListReader el(' ', 1);
        graph = el.read(clParser.getFileName().value());
    } else if (fileFormat == "mtx") {
        NetworKit::MatrixMarketReader mmr;
        auto res = mmr.read(clParser.getFileName().value());
        graph = matrixToUndirectedGraphSafe(res);
    } else {
        std::cerr << "Unknown file format (or wrong instance type specified)." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (clParser.getInstType() == "mc") {
        auto instance = sms::MaxCut(graph, clParser.getSeed());
        return instance;
    }

    if (clParser.getInstType() == "fi") {
        // ToDo: Seed not used!
        auto fi = sms::FrustrationIndex(graph);
        if (not fi.isConsistent()) {
            std::cerr << "Instance type is fi, but graph is not a signed graph" << std::endl;
            exit(EXIT_FAILURE);
        }
        return sms::frustrationIndexToMaxCut(fi);
    }

    // If we end up here, none of the above matched and we error out
    std::cerr << "Either wrong file format or instance type." << std::endl;
    exit(EXIT_FAILURE);
}

void solveInstance(const CLParser &clArgs, const sms::MaxCut &maxCutInstance) {
    nlohmann::ordered_json stats;

    // graph stats
    stats["graph stats"] = sms::graphStatsJson(maxCutInstance.getGraph());
    stats["graph stats"]["scaling factor"] = maxCutInstance.getScalingFactor();
    solveInstanceCore(maxCutInstance, clArgs, stats);

    if (clArgs.getStatsFile().has_value()) {
        auto statsFile = std::ofstream(clArgs.getStatsFile().value());
        statsFile << stats.dump(4) << std::endl;
    }
    std::cout << stats.dump(4) << std::endl;
}

/** main function for sms example */
int main(int argc, char **argv) {
    Aux::setNumberOfThreads(1); // Force NetworKit to only use one thread

    // start messages
    printCommandLineFlags(argc, argv, std::cout);
    printStartMessage();

    // parse args, read in problem and transform to MaxCut
    const auto clArgs = CLParser(argc, argv);
    auto maxCutInstance = readInstance(clArgs);
    sms::outputGraphStats(maxCutInstance.getGraph(), std::cout, true);
    maxCutInstance.scale();
    std::cout << "Scaling is " << maxCutInstance.getScalingFactor()
              << (maxCutInstance.allSolutionsHaveIntegerWeight() ? " and all solutions have integer values."
                                                                 : " and solutions might have fractional values.")
              << std::endl;
    solveInstance(clArgs, maxCutInstance);

    return EXIT_SUCCESS;
}
