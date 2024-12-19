sms_add_test(data_reducer_mc_gtest basic
        data_reducer_mc_gtest.cpp
        ../../instance/mc_solution.cpp
        ../../presol/data_reducer_mc.cpp
        ../../io/io.cpp
        ../../graph/gomory_hu_tree.cpp
        ../../graph/dinic.cpp
        ../../graph/fast_vertex_separator.cpp
        ../../graph/graphs.cpp
        ../../solver/enumeration_solver.cpp
        ../../solver/abstract_solver.cpp
        ../../solver/special_structures_solver.cpp
        ../../graph/bipartite.cpp
)

sms_add_test(presolver_mc_gtest basic
        presolver_mc_gtest.cpp
        ../../graph/biconnected_partition.cpp
        ../../graph/dinic.cpp
        ../../graph/fast_vertex_separator.cpp
        ../../graph/graphs.cpp
        ../../graph/gomory_hu_tree.cpp
        ../../presol/data_reducer_mc.cpp
        ../../io/io.cpp
        ../presolver_mc.cpp
        ../../solver/enumeration_solver.cpp
        ../../solver/special_structures_solver.cpp
        ../../solver/abstract_solver.cpp
        ../../graph/bipartite.cpp
)
