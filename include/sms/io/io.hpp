#ifndef SMS_IO_HPP
#define SMS_IO_HPP

#include <fstream>
#include <optional>
#include <regex>
#include <tuple>
#include <unistd.h>
#include <vector>

#include "networkit/graph/Graph.hpp"

#include "sms/instance/qubo.hpp"

NetworKit::Graph test_graph();

/*
 * Strict string to int parsing.
 * Errors out, if the string is not a valid integer.
 * Even if the string might have a valid prefix.
 */
int strictStringToInt(const std::string &input);

/*
 * Strict string to double parsing.
 * Errors out, if the string is not a valid double.
 * Even if the string might have a valid prefix.
 * Excepts C++ style double representation and 'D' for exponent.
 */
double strictStringToDouble(std::string input);

/* Temporary file with RAII
 * The file is closed when the object goes out of scope.
 */
class TmpFile {
public:
    TmpFile(const TmpFile &) = delete;
    TmpFile(TmpFile &&) = delete;
    TmpFile &operator=(const TmpFile &) = delete;
    TmpFile &operator=(TmpFile &&) = delete;

    explicit TmpFile(const std::string &nameTemplate = "/tmp/sms_tmpfile_XXXXXX");

    ~TmpFile() { close(fd_); }

    const std::string &getFilename() const { return name_; }

private:
    std::string name_;
    int fd_;
};

/**
 * Interprets the @a matrix as adjacency matrix of a graph. If @a matrix is non-symmetric, the graph
 * will be directed.
 * @param matrix
 * @return The graph having an adjacency matrix equal to @a matrix.
 */
template <class Matrix>
NetworKit::Graph matrixToUndirectedGraphSafe(const Matrix &matrix) {
    NetworKit::Graph G(std::max(matrix.numberOfRows(), matrix.numberOfColumns()), true, false);
    matrix.forNonZeroElementsInRowOrder([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        if (weight != 1)
            throw std::runtime_error("Weighted graphs are not supported here.");
        // assert(!G.hasEdge(u, v));
        if (!G.hasEdge(u, v))
            G.addEdge(u, v, weight);
    });

    return G;
}

class McObj {
public:
    unsigned int edges;
    unsigned int nodes;
    std::vector<std::string> comments;
    std::vector<std::tuple<int, int, double>> edge_list;

    McObj(unsigned int, unsigned int, const std::vector<std::tuple<int, int, double>> &,
          const std::vector<std::string> &);
};

McObj fileToMcObj(const std::string &path);

class BqObj {
public:
    unsigned int edges;
    unsigned int nodes;
    std::vector<std::string> comments;
    std::vector<std::tuple<int, int, double>> edge_list;

    BqObj(unsigned int, unsigned int, const std::vector<std::tuple<int, int, double>> &,
          const std::vector<std::string> &);
};

BqObj fileToBqObj(const std::string &path);

NetworKit::Graph bqFileToGraph(const std::string &path);

bool sortByBoth(const std::tuple<int, int, double> &a, const std::tuple<int, int, double> &b);

std::vector<std::string> fileToVector(const std::string &path);

std::tuple<std::vector<std::string>, std::vector<std::string>>
parseComments(const std::vector<std::string> &file_vector);

std::tuple<unsigned int, unsigned int, std::vector<std::string>> parseHeader(std::vector<std::string> file_vector);

std::vector<std::tuple<int, int, double>> parseEdgelist(const std::vector<std::string> &edgelist, int lower, int upper);

std::tuple<int, int, double> parseEdge(const std::string &edge_string);

bool checkMcEdgelist(const std::vector<std::tuple<int, int, double>> &edgelist);

bool checkBqEdgelist(const std::vector<std::tuple<int, int, double>> &edgelist);

std::vector<std::string> splitString(std::string input, char delimiter);

NetworKit::Graph edgelistToGraph(int nodes, const std::vector<std::tuple<int, int, double>> &edges);

NetworKit::Graph mcFileToGraph(const std::string &path);

#endif // SMS_IO_HPP
