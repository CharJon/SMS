sms_add_test(biconnected_partition_gtest basic biconnected_partition_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../graph/biconnected_partition.cpp
)

sms_add_test(bipartite_gtest basic bipartite_gtest.cpp
        ../../graph/bipartite.cpp
)

sms_add_test(gomory_hu_tree_gtest basic gomory_hu_tree_gtest.cpp
        ../gomory_hu_tree.cpp
        ../dinic.cpp
        ../../io/io.cpp
)


sms_add_test(graph_statistics_gtest basic graph_statistics_gtest.cpp
        ../../graph/graphs_statistics.cpp
        ../../graph/small_ccs.cpp
        ../../graph/graphs.cpp
        ../../graph/bipartite.cpp)

sms_add_test(dinic_gtest basic dinic_gtest.cpp
        ../dinic.cpp
)

sms_add_test(fast_vertex_separator_gtest basic fast_vertex_separator_gtest.cpp
        ../fast_vertex_separator.cpp
)

sms_add_test(graphs_gtest basic graphs_gtest.cpp
        ../../graph/graphs.cpp
        ../../io/io.cpp
)

sms_add_test(small_ccs_gtest basic small_ccs_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(odd_closed_walk_gtest scip odd_closed_walk_gtest.cpp
        ../../graph/odd_closed_walk.cpp
)

sms_add_test(ocw_separator_gtest basic ocw_separator_gtest.cpp
        ../ocw_seperator.cpp
        ../odd_closed_walk.cpp
)

sms_add_test(two_separator_gtest scip two_separator_gtest.cpp
        ../../graph/two_separator.cpp
)

sms_add_test(compact_graph_gtest scip compact_graph_gtest.cpp
)

sms_add_test(mst_heuristic_gtest scip mst_heuristic_gtest.cpp
        ../../graph/mst_heuristic.cpp
        ../../instance/mc_solution.cpp
        ../../graph/odd_closed_walk.cpp
        ../../auxiliary/union_find.cpp
        ../../graph/graphs.cpp
)

sms_add_test(gridded_ising_ccs_gtest scip gridded_ising_ccs_gtest.cpp
        ../../graph/gridded_ising_ccs.cpp
        ../../instance/ising.cpp
        ../../instance/gridded_ising.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp)
