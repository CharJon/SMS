sms_add_test(enumeration_solver_gtest basic enumeration_solver_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../instance/mc_solution.cpp
        ../../solver/enumeration_solver.cpp
        ../../solver/abstract_solver.cpp
        ../../graph/graphs.cpp
)

sms_add_test(enumeration_solver_gray_gtest basic enumeration_solver_gray_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../instance/mc_solution.cpp
        ../../solver/enumeration_solver_gray.cpp
        ../../solver/abstract_solver.cpp
        ../../graph/graphs.cpp
)

sms_add_test(gurobi_quadratic_solver_gtest gurobi
        gurobi_quadratic_solver_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../solver/gurobi_quadratic_solver.cpp
        ../../instance/qubo.cpp
        ../../instance/maxcut.cpp
        ../../instance/qubo_solution.cpp
        ../../io/qplib_parser.cpp
        ../../solver/abstract_solver.cpp
        ../../graph/graphs.cpp
)

sms_add_test(heuristic_solver_gtest basic
        heuristic_solver_gtest.cpp
        ../../graph/graphs.cpp
        ../../solver/abstract_solver.cpp
        ../../solver/heuristic_solver.cpp
        ../../auxiliary/mqlib.cpp
)

sms_add_test(special_structures_solver_gtets basic
        special_structures_solver_gtets.cpp
        ../../solver/special_structures_solver.cpp
        ../../graph/dinic.cpp
        ../../graph/graphs.cpp
        ../../graph/bipartite.cpp
        ../../instance/maxcut.cpp
        ../../solver/abstract_solver.cpp
        ../../graph/graphs.cpp
)
