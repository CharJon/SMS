sms_add_test(active_elements_gtest basic active_elements_gtest.cpp
)

sms_add_test(chrono_io_gtest basic chrono_io_gtest.cpp
)

sms_add_test(gray_code_gtest basic gray_code_gtest.cpp
)

sms_add_test(math_gtest basic math_gtest.cpp
)

sms_add_test(sms_statistics_gtest basic statistics_gtest.cpp
        ../../io/io.cpp
        ../../graph/graphs.cpp
        ../../graph/small_ccs.cpp
)

sms_add_test(subset_enumerator_gtest basic subset_enumerator_gtest.cpp
        ../../auxiliary/subset_enumerator.cpp
)

sms_add_test(union_find_gtest basic union_find_gtest.cpp
        ../../auxiliary/union_find.cpp
)
