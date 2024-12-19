add_executable(bench_cache cache.cpp)

target_include_directories(bench_cache
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        )

target_link_libraries(bench_cache
        benchmark::benchmark)
