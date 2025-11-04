#include "test_common.hpp"

namespace fp {
namespace test {

void run_divide_tests() {
    using q16 = q<1, 15, fp::test::Backend>;
    using q8  = q<1, 7, fp::test::Backend>;
    using q32 = q<3, 29, fp::test::Backend>;
    using q16_ref = q<1, 15, fp::test::Backend>;
    using q32_ref = q<3, 29, fp::test::Backend>;

    std::puts("\n--- Division Tests ---");

    // Basic divide tests with matching operand widths
    // Q1.15 can represent [-1.0, 1.0), so 0.375/0.50 = 0.75 fits
    check_div<q16,  q16,  1,15>("16÷16->16 (fp::test::Backend prio-1)", 0.375f, 0.50f);
    // Q1.7 can represent [-1.0, 1.0), so 0.25/0.50 = 0.5 fits
    check_div<q8,   q8,   1,7> ("8÷8->8  (fp::test::Backend prio-2)",  0.25f, 0.50f);
    // Q3.29 can represent up to ~4.0, so 0.75/0.50 = 1.5 fits
    check_div<q16,  q16,  3,29>("16÷16->32 (fallback)",            0.75f, 0.50f);

    // Mixed width operands (results fit in Q1.15 and Q1.7)
    check_div<q16,  q8,   1,15>("16÷8->16",  0.30f, 0.60f);  // 0.5
    check_div<q8,   q16,  1,7> ("8÷16->8",   0.20f, 0.40f);  // 0.5

    // Negative values (results fit in Q1.15 range)
    check_div<q16,  q16,  1,15>("Signed: (-0.375)/0.5",  -0.375f, 0.50f);  // -0.75
    check_div<q16,  q16,  1,15>("Signed: 0.375/(-0.5)",   0.375f, -0.50f);  // -0.75
    check_div<q16,  q16,  1,15>("Signed: (-0.375)/(-0.5)", -0.375f, -0.50f);  // 0.75

    // Test operator/ with values that fit
    {
        auto c = q16::from_float(0.375f);
        auto d = q16::from_float(0.50f);
        auto quotient = c / d;  // Should use operator/
        float got = quotient.to_float();
        float want = 0.375f / 0.50f;  // = 0.75
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Operator / test", got, want, 2.0f * lsb);
    }

    // Test free function div_as with result that needs more integer bits
    {
        auto e = q16::from_float(0.75f);
        auto f = q16::from_float(0.50f);
        auto quotient2 = fp::div_as<3,29>(e, f);  // Output Q3.29 (can hold 1.5)
        float got = quotient2.to_float();
        float want = 0.75f / 0.50f;  // = 1.5
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        expect_near("div_as<3,29> test", got, want, 2.0f * lsb);
    }

    // Test Reference backend
    check_div<q16_ref, q16_ref, 1,15>("Reference 16÷16->16", 0.375f, 0.50f);
    {
        auto a = q32_ref::from_float(2.0f);
        auto b = q32_ref::from_float(1.0f);
        auto c = a / b;
        float got = c.to_float();
        float want = 2.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        expect_near("Reference 32÷32->32", got, want, 2.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_divide_tests();
    return 0;
}
#endif
