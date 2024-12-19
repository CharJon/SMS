cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 20)

project(sms)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native")

set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "") # Enable/Disable output of compile commands during generation

# New linker is broken on apple with gcc
if (APPLE)
    add_link_options("-Wl,-ld_classic")
endif ()

# find GUROBI is desired
option(GUROBI "Build Gurobi part of the project" OFF)

if (GUROBI)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
    find_package(GUROBI REQUIRED)
    include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
endif ()

# BLAS
if (APPLE)
    #set(BLA_VENDOR "OpenBLAS") <- this should? be the way, but does not work
    # we hardcode homebrew
    find_package(OpenBLAS REQUIRED PATHS "/opt/homebrew/opt/openblas")
    set(BLAS_LIBRARIES "${OpenBLAS_LIBRARIES}")
    include_directories("${OpenBLAS_INCLUDE_DIRS}")
else ()
    find_package(BLAS REQUIRED)
endif ()


### This Project #####
add_subdirectory(src)

# NO hpp files here!
set(SOURCES
        src/auxiliary/bipartition.cpp
        src/auxiliary/chrono_io.cpp
        src/auxiliary/cl_parser.cpp
        src/auxiliary/mqlib.cpp
        src/auxiliary/scip.cpp
        src/auxiliary/union_find.cpp
        src/auxiliary/subset_enumerator.cpp
        src/branch/branchrule_degree.cpp
        src/branch/branchrule_nm.cpp
        src/conshdlr/cliques.cpp
        src/conshdlr/conshdlr_triangles.cpp
        src/conshdlr/conshdlr_holes.cpp
        src/conshdlr/conshdlr_cliques.cpp
        src/conshdlr/conshdlr_oscw.cpp
        src/disp/disp.cpp
        src/eventhdlr/eventhdlr_history.cpp
        src/eventhdlr/eventhdlr_root.cpp
        src/graph/biconnected_partition.cpp
        src/graph/bipartite.cpp
        src/graph/dinic.cpp
        src/graph/fast_vertex_separator.cpp
        src/graph/gomory_hu_tree.cpp
        src/graph/graphs.cpp
        src/graph/small_ccs.cpp
        src/graph/two_separator.cpp
        src/graph/graphs_statistics.cpp
        src/graph/mst_heuristic.cpp
        src/graph/odd_closed_walk.cpp
        src/graph/ocw_seperator.cpp
        src/instance/convert.cpp
        src/instance/maxcut.cpp
        src/instance/mc_solution.cpp
        src/instance/ising.cpp
        src/instance/gridded_ising.cpp
        src/instance/qubo.cpp
        src/instance/frustration_index.cpp
        src/io/io.cpp
        src/io/dimacs_reader.cpp
        src/io/ising_reader.cpp
        src/io/json.cpp
        src/io/mc_read_write.cpp
        src/io/qplib_parser.cpp
        src/pheur/kl_heuristic.cpp
        src/pheur/pheur_kl.cpp
        src/pheur/pheur_mst.cpp
        src/pheur/pheur_mqlib.cpp
        src/presol/data_reducer_mc.cpp
        src/presol/presolver_mc.cpp
        src/probdata/mc_rooted.cpp
        src/probdata/mc_edges_only.cpp
        src/solver/basic_ocw.cpp
        src/solver/abstract_solver.cpp
        src/solver/enumeration_solver.cpp
        src/solver/heuristic_solver.cpp
        src/solver/special_structures_solver.cpp
        src/conshdlr/maximal_cliques.cpp
        src/pheur/burer_heuristic.cpp
        src/pheur/pheur_burer.cpp
)
if (GUROBI)
    list(APPEND SOURCES src/solver/gurobi_quadratic_solver.cpp)
endif ()

add_executable(sms
        src/mains/main.cpp
        ${SOURCES}
)

target_include_directories(sms
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
if (GUROBI)
    target_include_directories(sms PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_compile_definitions(sms PRIVATE SMS_GUROBI)
endif ()

### Checks exe ####

set(SOURCES_CHECKS
        src/mains/checks.cpp
        src/graph/graphs.cpp
        src/graph/small_ccs.cpp
        src/io/io.cpp
        src/io/ising_reader.cpp
        src/io/dimacs_reader.cpp
        src/io/qplib_parser.cpp
        src/instance/ising.cpp
        src/instance/gridded_ising.cpp
        src/instance/maxcut.cpp
        src/instance/qubo.cpp)

# An executable for checking the correctness of input files
add_executable(checks ${SOURCES_CHECKS})

target_include_directories(checks
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)


### Extern Projects #####

# NetworKIT
if (NOT NETWORKIT_DIR)
    message("No NETWORKIT_DIR set via CMAKE, will try environment variable NETWORKIT_DIR: $ENV{NETWORKIT_DIR}")
    set(NETWORKIT_DIR "$ENV{NETWORKIT_DIR}" CACHE PATH "Path to NetworKIT directory" FORCE)
endif ()
if (NOT NETWORKIT_DIR)
    add_subdirectory(extern/networkit)
else ()
    add_library(networkit SHARED IMPORTED)
    set_target_properties(networkit PROPERTIES
            IMPORTED_LOCATION ${NETWORKIT_DIR}/lib/libnetworkit.so
            INTERFACE_INCLUDE_DIRECTORIES ${NETWORKIT_DIR}/include
    )
endif ()

# SCIP
if (NOT SCIP_DIR)
    message("No SCIP_DIR set via CMAKE, will try environment variable SCIP_DIR: $ENV{SCIP_DIR}")
    set(SCIP_DIR "$ENV{SCIP_DIR}" CACHE PATH "Path to SCIP directory" FORCE)
else ()
    message(STATUS "Using SCIP library from here: ${SCIP_DIR}.")
endif ()
find_package(SCIP PATHS ${SCIP_DIR} NO_DEFAULT_PATH REQUIRED)

# json
add_subdirectory(extern/json)

#cxxopts
add_subdirectory(extern/cxxopts)

# mqlib
add_subdirectory(extern/MQLib)

# gtest
option(TESTS "Build the tests" ON)
if (TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/gtest" "extern/googletest")

    # This macro has the testname as the first parameter and the type as the second.
    # All further arguments have to be .cpp files containing the relevant sources
    macro(sms_add_test TESTNAME TEST_TYPE)
        if ((NOT GUROBI) AND ("${TEST_TYPE}" STREQUAL "gurobi"))
            message(WARNING "Test ignored, is flagged gurobi ${TEST_TYPE}.")
        else ()
            # create an executable in which the tests will be stored
            add_executable(${TESTNAME} ${ARGN})
            # link the Google test infrastructure, mocking library, and a default main function to
            # the test executable.  Remove g_test_main if writing your own main function.
            target_include_directories(${TESTNAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

            if ("${TEST_TYPE}" STREQUAL "basic")
                target_link_libraries(${TESTNAME} gtest gmock gtest_main ${BLAS_LIBRARIES} networkit mqlib nlohmann_json::nlohmann_json cxxopts)
            elseif ("${TEST_TYPE}" STREQUAL "scip")
                target_link_libraries(${TESTNAME} gtest gmock gtest_main ${BLAS_LIBRARIES} networkit mqlib nlohmann_json::nlohmann_json cxxopts ${SCIP_LIBRARIES})
            elseif ("${TEST_TYPE}" STREQUAL "gurobi")
                target_link_libraries(${TESTNAME} gtest gmock gtest_main ${BLAS_LIBRARIES} networkit mqlib nlohmann_json::nlohmann_json cxxopts ${SCIP_LIBRARIES})
                target_include_directories(${TESTNAME} PRIVATE ${GUROBI_INCLUDE_DIRS})
                target_link_libraries(${TESTNAME}
                        optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY}
                        ${GUROBI_LIBRARY})
            else ()
                message(SEND_ERROR "Test type unknown, got ${TEST_TYPE}.")
            endif ()

            # gtest_discover_tests replaces gtest_add_tests,
            # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
            gtest_discover_tests(${TESTNAME}
                    # set a working directory so your project root so that you can find test data via paths relative to the project root
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
            )
            set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
        endif ()
    endmacro()

    # ADD TEST DIRECTORIES HERE
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/auxiliary/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/branch/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/conshdlr/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/graph/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/pheur/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/presol/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/probdata/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/io/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/instance/test)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/solver/test)
endif ()

# gbenchmmark
option(MICROBENCH "Build microbenchmarks" OFF)
if (MICROBENCH)
    set(BENCHMARK_ENABLE_TESTING OFF)
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/gbenchmark" "extern/gbenchmark")
    add_subdirectory(micro_benchmark)
endif ()


## This project links
target_link_libraries(sms
        ${SCIP_LIBRARIES}
        ${BLAS_LIBRARIES}
        cxxopts
        networkit
        nlohmann_json::nlohmann_json
        mqlib
)
if (GUROBI)
    target_link_libraries(sms
            optimized ${GUROBI_CXX_LIBRARY}
            debug ${GUROBI_CXX_DEBUG_LIBRARY}
            ${GUROBI_LIBRARY}
    )
endif ()

target_link_libraries(checks
        ${BLAS_LIBRARIES}
        cxxopts
        networkit
        nlohmann_json::nlohmann_json
)


## install
message(STATUS ${CMAKE_INSTALL_PREFIX})
install(TARGETS sms RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})


# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()



