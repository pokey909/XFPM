#include "test_common.hpp"

namespace fp {
namespace test {

void run_multiply_tests() {
    using q16 = q<1, 15, fp::test::Backend>;
    using q8  = q<1, 7, fp::test::Backend>;
    using q32 = q<3, 29, fp::test::Backend>;
    using q16_ref = q<1, 15, fp::test::Backend>;
    using q32_ref = q<3, 29, fp::test::Backend>;

    std::puts("\n--- Multiply Tests ---");

    // Basic multiply tests with matching operand widths
    check_mul<q8,  q8,  1,7> ("8x8->8  (fp::test::Backend prio-2)",  0.5f, 0.75f);
    check_mul<q8,  q8,  1,15>("8x8->16 (fallback)",               0.5f, 0.75f);
    check_mul<q16, q16, 1,15>("16x16->16 (fp::test::Backend prio-1)", 0.5f, 0.75f);
    check_mul<q16, q16, 3,29>("16x16->32 (fallback)",             0.5f, 0.75f);

    // Arbitrary operand widths (bucket logic)
    check_mul<q<6,6,fp::test::Backend>, q<6,6,fp::test::Backend>, 1,15>("12b×12b -> 16b bucket", 0.15f, 0.80f);

    // Mixed width operands (results fit in Q1.15)
    check_mul<q8,  q16, 1,15>("8b×16b -> 16b", 0.30f, 0.80f);
    check_mul<q16, q8,  1,15>("16b×8b -> 16b", 0.30f, 0.80f);
    check_mul<q16, q32, 3,29>("16b×32b -> 32b", 0.30f, 0.80f);

    // Explicit result format via mul_as (shift the fractional bits)
    {
        auto a = q16::from_float(0.5f);
        auto b = q16::from_float(0.5f);
        auto c = fp::mul_as<3,10>(a, b);  // Q3.10 instead of Q1.15
        float got = c.to_float();
        float want = 0.5f * 0.5f;  // = 0.25
        const float lsb = 1.0f / static_cast<float>(1u << 10);
        expect_near("16x16 -> Q3.10 (shift test)", got, want, 2.0f * lsb);
    }

    // Test Reference backend
    check_mul<q16_ref, q16_ref, 1,15>("Reference 16x16->16", 0.5f, 0.75f);
    check_mul<q16_ref, q16_ref, 3,29>("Reference 16x16->32", 0.5f, 0.75f);

    // Test edge case: near saturation
    {
        auto a = q16::from_float(0.999f);
        auto b = q16::from_float(1.0f);
        auto c = a * b;
        float got = c.to_float();
        float want = 0.999f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Saturation-safe multiply near max (Q1.15)", got, want, 2.0f * lsb);
    }

    // Test different Q formats
    {
        using q3_13 = q<3, 13, fp::test::Backend>;
        auto a = q3_13::from_float(2.0f);
        auto b = q16::from_float(0.2f);
        auto c = fp::mul_as<3,13>(a, b);
        float got = c.to_float();
        float want = 2.0f * 0.2f;  // = 0.4
        const float lsb = 1.0f / static_cast<float>(1u << 13);
        expect_near("Q3.13 * Q1.15 -> Q3.13", got, want, 2.0f * lsb);
    }

    // Test signed values
    {
        auto a = q16::from_float(-0.5f);
        auto b = q16::from_float(0.75f);
        auto c = a * b;
        float got = c.to_float();
        float want = -0.5f * 0.75f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Signed: (-0.5)*0.75", got, want, 2.0f * lsb);
    }

    {
        auto a = q16::from_float(-0.3f);
        auto b = q16::from_float(0.6f);
        auto c = a * b;
        float got = c.to_float();
        float want = -0.3f * 0.6f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Signed: (-0.3)*0.6", got, want, 3.0f * lsb);
    }

    // Simple test to print
    {
        auto x = q16::from_float(0.5f);
        auto y = q16::from_float(0.75f);
        auto z = x * y;
        std::printf("%f  %f\n", z.to_float(), 0.5f * 0.75f);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_multiply_tests();
    return 0;
}
#endif
