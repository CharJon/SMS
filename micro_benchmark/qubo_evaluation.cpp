#include <iostream>
#include <random>
#include <vector>

#include "benchmark/benchmark.h"

#include "sms/auxiliary/math.hpp"

static void quboEvalBench(benchmark::State &s) {
    int dim = s.range(0);

    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()}; // Generates random integers

    std::uniform_int_distribution<int> dist{0, 1};
    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };

    std::vector<double> x(dim);
    std::generate(begin(x), end(x), gen);

    std::vector<double> Q(dim * dim);
    std::uniform_real_distribution<double> dist_q{-10000.0, 10000.0};
    auto gen_indices = [&dist_q, &mersenne_engine]() { return dist_q(mersenne_engine); };
    std::generate(begin(Q), end(Q), gen_indices);

    int iter = static_cast<int>(pow(2.0, static_cast<double>(dim - 1)));

    for (auto _ : s) {
        for (int i = 0; i < iter; i++)
            benchmark::DoNotOptimize(sms::quboEvaluation(&Q[0], &x[0], dim, 1.0));
    }

    s.counters["Variables"] = dim;
    s.counters["Runs"] = iter;
    s.counters["Runtime"] =
        benchmark::Counter(iter, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kAvgThreadsRate
                                     | benchmark::Counter::kInvert);
}

// Register the function as a benchmark
BENCHMARK(quboEvalBench)->ArgsProduct({benchmark::CreateDenseRange(10, 25, 1)});

// Run the benchmark
BENCHMARK_MAIN();
