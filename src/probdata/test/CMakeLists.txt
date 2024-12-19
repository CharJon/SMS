sms_add_test(prob_data_mc_gtest scip
        prob_data_mc_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../auxiliary/scip.cpp
        ../../instance/maxcut.cpp
        ../../probdata/prop_data_mc.cpp
)

sms_add_test(probdata_mc_rooted_gtest scip probdata_mc_rooted_gtest.cpp
        ../../auxiliary/mqlib.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/graphs_statistics.cpp
        ../../graph/small_ccs.cpp
        ../../instance/mc_solution.cpp
        ../../auxiliary/scip.cpp
        ../../auxiliary/union_find.cpp
        ../../graph/bipartite.cpp
        ../../graph/mst_heuristic.cpp
        ../../graph/odd_closed_walk.cpp
        ../../graph/ocw_seperator.cpp
        ../../branch/branchrule_degree.cpp
        ../../conshdlr/conshdlr_oscw.cpp
        ../../conshdlr/conshdlr_triangles.cpp
        ../../conshdlr/conshdlr_holes.cpp
        ../../eventhdlr/eventhdlr_history.cpp
        ../../eventhdlr/eventhdlr_root.cpp
        ../../pheur/pheur_mst.cpp
        ../../pheur/pheur_mqlib.cpp
        ../../probdata/mc_rooted.cpp
        ../../solver/heuristic_solver.cpp
        ../../solver/abstract_solver.cpp
)
