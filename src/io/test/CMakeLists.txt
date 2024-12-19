sms_add_test(io_gtest scip io_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
        ../../instance/qubo.cpp
)

sms_add_test(dimacs_reader_gtest scip dimacs_reader_gtest.cpp
        ../../io/io.cpp
        ../../io/dimacs_reader.cpp
        ../../instance/ising.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(ising_reader_gtest scip ising_reader_gtest.cpp
        ../../io/ising_reader.cpp
        ../../io/io.cpp
        ../../instance/ising.cpp
        ../../instance/gridded_ising.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(mc_read_write_gtest basic mc_read_write_gtest.cpp
        ../../io/mc_read_write.cpp
        ../../io/io.cpp
)

sms_add_test(qplib_parser_gtest scip qplib_parser_gtest.cpp
        ../../io/qplib_parser.cpp
        ../../io/io.cpp
        ../../instance/qubo.cpp
)