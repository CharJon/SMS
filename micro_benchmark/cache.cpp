#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "benchmark/benchmark.h"

static void cacheBench(benchmark::State &s) {
    int bytes = 1 << s.range(0);
    int count = (bytes / sizeof(int)) / 2;
    bool random = s.range(1);

    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine{rnd_device()}; // Generates random integers

    std::uniform_int_distribution<int> dist{std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
    auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };

    std::vector<int> v(count);
    std::generate(begin(v), end(v), gen);

    std::vector<int> indices(count);
    if (random) {
        std::uniform_int_distribution<int> dist_indices{0, count};
        auto gen_indices = [&dist_indices, &mersenne_engine]() { return dist_indices(mersenne_engine); };
        std::generate(begin(indices), end(indices), gen_indices);
    } else {
        for (int i = 0; i < count; i++) {
            indices[i] = i;
        }
    }

    for (auto _ : s) {
        uint64_t sum = 0;
        for (int i : indices)
            sum += v[i];
        benchmark::DoNotOptimize(sum);
    }

    s.SetBytesProcessed(long(s.iterations()) * long(bytes));
    s.SetLabel(std::to_string(bytes / 1024) + "kB");
}

// Register the function as a benchmark
BENCHMARK(cacheBench)->ArgsProduct({benchmark::CreateDenseRange(17, 25, 1), {0, 1}});

// Function to be benchmarked
double sumSelected(const std::vector<uint8_t > &bools, const std::vector<double> &values) {
    double sum = 0.0;
    for (size_t i = 0; i < bools.size(); ++i) {
        if (bools[i]) {
            sum += values[i];
        } else {
            sum -= values[i];
        }
    }
    return sum;
}

// Function to be benchmarked
double sumSelected2(const std::vector<uint8_t > &bools, const std::vector<double> &values) {
    double sum = 0.0;
    for (size_t i = 0; i < bools.size(); ++i) {
        sum += (1 - 2 * bools[i]) * values[i];
    }
    return sum;
}

// Function to be benchmarked
double sumSelected3(const std::vector<uint8_t > &bools, const std::vector<double> &values) {
    double sum = 0.0;
    for (size_t i = 0; i < bools.size(); ++i) {
        sum += bools[i] ? values[i] : -values[i];
    }
    return sum;
}

// Benchmark for the sumSelected function
static void BM_SumSelected(benchmark::State &state) {
    // Generate random vectors
    size_t vector_size = state.range(0); // Control the size via benchmark range
    std::vector<uint8_t > bools(vector_size);
    std::vector<double> values(vector_size);

    // Random number generator
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.0, 100.0);

    // Fill vectors with random values
    for (size_t i = 0; i < vector_size; ++i) {
        bools[i] = gen() % 2; // Set some pattern for bools (true/false)
        values[i] = dis(gen); // Random doubles
    }

    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(sumSelected(bools, values));
    }
}

// Benchmark for the sumSelected function
static void BM_SumSelected2(benchmark::State &state) {
    // Generate random vectors
    size_t vector_size = state.range(0); // Control the size via benchmark range
    std::vector<uint8_t > bools(vector_size);
    std::vector<double> values(vector_size);

    // Random number generator
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.0, 100.0);

    // Fill vectors with random values
    for (size_t i = 0; i < vector_size; ++i) {
        bools[i] = gen() % 2; // Set some pattern for bools (true/false)
        values[i] = dis(gen); // Random doubles
    }

    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(sumSelected2(bools, values));
    }
}

// Benchmark for the sumSelected function
static void BM_SumSelected3(benchmark::State &state) {
    // Generate random vectors
    size_t vector_size = state.range(0); // Control the size via benchmark range
    std::vector<uint8_t > bools(vector_size);
    std::vector<double> values(vector_size);

    // Random number generator
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(0.0, 100.0);

    // Fill vectors with random values
    for (size_t i = 0; i < vector_size; ++i) {
        bools[i] = gen() % 2; // Set some pattern for bools (true/false)
        values[i] = dis(gen); // Random doubles
    }

    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(sumSelected3(bools, values));
    }
}

// Register the benchmark with different input sizes
BENCHMARK(BM_SumSelected)->Range(1 << 10, 1 << 20);  // From 1024 to 1048576 elements
BENCHMARK(BM_SumSelected2)->Range(1 << 10, 1 << 20); // From 1024 to 1048576 elements
BENCHMARK(BM_SumSelected3)->Range(1 << 10, 1 << 20); // From 1024 to 1048576 elements

// Run the benchmark
BENCHMARK_MAIN();
