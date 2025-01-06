#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "benchmark/benchmark.h"

template <typename T>
struct Edge
{
    T u;
    T v;
    double weight;
};

template <typename T>
static void cacheBench(benchmark::State& s)
{
    bool random = s.range(1);

    int numElements = 1 << s.range(0); //bytes / sizeof(Edge<T>);
    int64_t bytes = numElements * sizeof(Edge<T>);

    std::vector<T> ids(numElements);
    std::iota(ids.begin(), ids.end(), 0); // Fill with 0, 1, ..., numElements-1
    if (random)
    {
        std::mt19937 mersenne_engine(100); // Generates random integers
        std::shuffle(ids.begin(), ids.end(), mersenne_engine);
    }

    std::vector<Edge<T>> edges(numElements);
    for (T i = 0; i < numElements - 1; ++i)
    {
        auto from = ids[i];
        auto to = ids[i + 1];
        edges[from] = Edge<T>(from, to, i);
    }
    auto last = ids.back();
    edges[last] = Edge<T>(last, -1, 0);

    for (auto _ : s)
    {
        uint64_t sum = 0;
        auto next = ids.front();
        while (next != -1)
        {
            sum += edges[next].weight;
            next = edges[next].v;
        }
        benchmark::DoNotOptimize(sum);
    }

    //s.SetBytesProcessed(static_cast<long>(s.iterations()) * static_cast<long>(bytes));
    s.SetLabel(std::to_string(bytes / 1024) + "kB");
}

// Register the function as a benchmark
BENCHMARK(cacheBench<int32_t>)->ArgsProduct({benchmark::CreateDenseRange(10, 20, 1), {0}});
BENCHMARK(cacheBench<int64_t>)->ArgsProduct({benchmark::CreateDenseRange(10, 20, 1), {0}});
BENCHMARK(cacheBench<int32_t>)->ArgsProduct({benchmark::CreateDenseRange(10, 20, 1), {1}});
BENCHMARK(cacheBench<int64_t>)->ArgsProduct({benchmark::CreateDenseRange(10, 20, 1), {1}});

// Run the benchmark
BENCHMARK_MAIN();
