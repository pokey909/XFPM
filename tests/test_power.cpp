#include "test_common.hpp"

namespace fp {
namespace test {

void run_power_tests() {
    using q16 = q<1, 15, fp::test::Backend>;
    using q8  = q<1, 7, fp::test::Backend>;

    std::puts("\n--- Power Tests ---");

    // Test case 1: pow(2.0, 3.0) = 8.0
    // Need Q format that can hold 8.0 for result, and 3.0 for exponent
    // Use Q4.12 for base (can hold up to ~16) and Q3.13 for exponent (can hold up to ~4)
    {
        using q4_12 = q<4,12, fp::test::Backend>;
        using q3_13 = q<3,13, fp::test::Backend>;  // Can hold [-4.0, 4.0)
        check_pow<q4_12, q3_13>("pow(2.0, 3.0) = 8.0", 2.0f, 3.0f, 2.0f);
    }

    // Test case 2: pow(0.5, 2.0) = 0.25
    // Use Q3.13 for exponent to hold 2.0
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q16, q3_13>("pow(0.5, 2.0) = 0.25", 0.5f, 2.0f, 2.0f);
    }

    // Test case 3: pow(4.0, 0.5) = 2.0 (square root of 4)
    {
        using q4_12 = q<4,12, fp::test::Backend>;
        check_pow<q4_12, q16>("pow(4.0, 0.5) = 2.0", 4.0f, 0.5f, 2.0f);
    }

    // Test case 4: pow(0.25, 2.0) = 0.0625
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q16, q3_13>("pow(0.25, 2.0) = 0.0625", 0.25f, 2.0f, 2.0f);
    }

    // Test case 5: pow(3.0, 2.0) = 9.0
    // Result 9.0 needs Q5.11 (can hold [-16.0, 16.0))
    {
        using q5_11 = q<5,11, fp::test::Backend>;
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q5_11, q3_13>("pow(3.0, 2.0) = 9.0", 3.0f, 2.0f, 2.0f);
    }

    // Test case 6: pow(0.5, 0.5) ≈ 0.707
    check_pow<q16, q16>("pow(0.5, 0.5) ≈ 0.707", 0.5f, 0.5f, 2.0f);

    // Test case 7: Fractional exponent - pow(8.0, 1/3) ≈ 2.0 (cube root)
    {
        using q4_12 = q<4,12, fp::test::Backend>;
        float one_third = 1.0f / 3.0f;
        check_pow<q4_12, q16>("pow(8.0, 1/3) ≈ 2.0", 8.0f, one_third, 2.0f);
    }

    // Test case 8: pow(x, 0) = 1 for any x
    // Exponent 0.0 fits in Q1.15
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q3_13, q16>("pow(2.0, 0.0) = 1.0", 2.0f, 0.0f, 2.0f);
    }

    // Test case 9: pow(1.0, x) = 1.0 for any x (but x must fit)
    // Exponent 3.0 needs Q3.13
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q16, q3_13>("pow(1.0, 3.0) = 1.0", 1.0f, 3.0f, 2.0f);
    }

    // Test case 10: Negative exponent - pow(2.0, -1.0) = 0.5
    // Base 2.0 needs Q3.13, exponent -1.0 fits in Q1.15
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        check_pow<q3_13, q16>("pow(2.0, -1.0) = 0.5", 2.0f, -1.0f, 2.0f);
    }

    // Test zero/negative base handling
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        auto base_zero = q16::from_float(0.0f);
        auto exp = q3_13::from_float(2.0f);  // Use Q3.13 for exponent 2.0
        auto result = base_zero.pow(exp);
        float got = result.to_float();
        expect_near("pow(0.0, 2.0) = 0.0 (error handling)", got, 0.0f, 0.001f);
    }

    {
        using q3_13 = q<3,13, fp::test::Backend>;
        auto base_neg = q16::from_float(-0.5f);
        auto exp = q3_13::from_float(2.0f);  // Use Q3.13 for exponent 2.0
        auto result = base_neg.pow(exp);
        float got = result.to_float();
        expect_near("pow(-0.5, 2.0) = 0.0 (error handling)", got, 0.0f, 0.001f);
    }

    // Test with 32-bit base (priority 1 path)
    // Result 8.0 needs Q5.27 (can hold [-16.0, 16.0))
    {
        using q5_27 = q<5,27, fp::test::Backend>;
        using q3_13 = q<3,13, fp::test::Backend>;
        auto base = q5_27::from_float(2.0f);
        auto exp = q3_13::from_float(3.0f);  // Use Q3.13 for exponent 3.0
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 27);
        expect_near("pow(2.0, 3.0) with Q32 base (Xtensa prio-1)", got, want, 2.0f * lsb);
    }

    // Test with 8-bit base (fallback to priority 0)
    {
        using q3_13 = q<3,13, fp::test::Backend>;
        auto base = q8::from_float(0.5f);
        auto exp = q3_13::from_float(2.0f);
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 0.25f;
        const float lsb = 1.0f / static_cast<float>(1u << 7);
        expect_near("pow(0.5, 2.0) with Q8 (fallback)", got, want, 2.0f * lsb);
    }

    // Test Reference backend parity
    // Base 2.0 and result 8.0 need Q4.12 (can hold [-8.0, 8.0))
    {
        using q4_12_ref = q<4,12, fp::test::Backend>;
        using q3_13_ref = q<3,13, fp::test::Backend>;
        auto base = q4_12_ref::from_float(2.0f);
        auto exp = q3_13_ref::from_float(3.0f);
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 12);
        expect_near("pow(2.0, 3.0) Reference backend", got, want, 10.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_power_tests();
    return 0;
}
#endif
