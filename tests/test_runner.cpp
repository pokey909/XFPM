#include "test_common.hpp"
#include <cstdio>

// Forward declarations of test functions
namespace fp {
namespace test {
    void run_multiply_tests();
    void run_divide_tests();
    void run_logarithm_tests();
    void run_antilogarithm_tests();
    void run_sqrt_tests();
    void run_power_tests();
    void run_array_ops_tests();
    void run_vector_ops_tests();
}
}

int main() {
    std::puts("===============================================");
    std::puts("  Fixed-Point Library Test Suite");
    std::puts("===============================================");

    // Run all test suites
    fp::test::run_multiply_tests();
    fp::test::run_divide_tests();
    fp::test::run_logarithm_tests();
    fp::test::run_antilogarithm_tests();
    fp::test::run_sqrt_tests();
    fp::test::run_power_tests();
    fp::test::run_array_ops_tests();
    fp::test::run_vector_ops_tests();

    // Summary
    std::puts("\n===============================================");
    if (fp::test::failures == 0) {
        std::puts("  ALL TESTS PASSED");
        std::puts("===============================================");
        return 0;
    } else {
        std::printf("  %d TEST(S) FAILED\n", fp::test::failures);
        std::puts("===============================================");
        return 1;
    }
}
