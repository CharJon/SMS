sms_add_test(instance_qubo_gtest basic
        qubo_gtest.cpp
        ../qubo.cpp
        ../../io/io.cpp
)

sms_add_test(instance_maxcut_gtest basic maxcut_gtest.cpp
        ../maxcut.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(convert_gtest basic convert_gtest.cpp
        ../../instance/convert.cpp
        ../../io/io.cpp
        ../../io/ising_reader.cpp
        ../../instance/qubo.cpp
        ../../instance/maxcut.cpp
        ../../instance/ising.cpp
        ../../instance/gridded_ising.cpp
        ../../instance/frustration_index.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../graph/gridded_ising_ccs.cpp
        ../../solver/enumeration_solver.cpp
        ../../solver/abstract_solver.cpp
)

sms_add_test(instance_ising_gtest basic ising_gtest.cpp
        ../ising.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(gridded_ising_gtest basic gridded_ising_gtest.cpp
        ../ising.cpp
        ../gridded_ising.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(frustration_index_graph_gtest basic frustration_index_gtest.cpp
        ../frustration_index.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp)

sms_add_test(mc_solution_gtest basic
        mc_solution_gtest.cpp
        ../../instance/mc_solution.cpp
)

sms_add_test(qubo_solution_gtest scip qubo_solution_gtest.cpp
        ../../instance/qubo_solution.cpp
        ../../instance/qubo.cpp
        ../../io/io.cpp
)
