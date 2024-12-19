sms_add_test(triangles_gtest scip triangles_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../graph/odd_closed_walk.cpp
)

sms_add_test(holes_gtest scip holes_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../graph/odd_closed_walk.cpp
)

sms_add_test(cliques_gtest scip cliques_gtest.cpp
        ../cliques.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../pheur/kl_heuristic.cpp
        ../../solver/enumeration_solver.cpp
        ../../solver/abstract_solver.cpp
        ../../auxiliary/bipartition.cpp
        ../../instance/mc_solution.cpp
        ../../auxiliary/subset_enumerator.cpp
)