#include "test_common.hpp"

namespace fp {
namespace test {

void run_sqrt_tests() {
    using q16 = q<1, 15, fp::test::Backend>;
    using q8  = q<1, 7, fp::test::Backend>;
    using q32 = q<3, 29, fp::test::Backend>;
    using q16_ref = q<1, 15, fp::test::Backend>;

    std::puts("\n--- Square Root Tests ---");

    // Basic sqrt tests with common values
    {
        auto x = q16::from_float(0.25f);  // sqrt(0.25) = 0.5
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sqrt(0.25)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16::from_float(0.5f);  // sqrt(0.5) ≈ 0.707
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sqrt(0.5)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16::from_float(0.16f);  // sqrt(0.16) = 0.4
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.16f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sqrt(0.16)", got, want, 2.0f * lsb);
    }

    // Test with 32-bit input
    {
        auto x = q32::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        expect_near("sqrt(0.25) with Q32", got, want, 2.0f * lsb);
    }

    {
        auto x = q32::from_float(2.0f);  // sqrt(2) ≈ 1.414
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(2.0f);
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        expect_near("sqrt(2.0) with Q32", got, want, 2.0f * lsb);
    }

    // Test with 8-bit input
    {
        auto x = q8::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 7);
        expect_near("sqrt(0.25) with Q8", got, want, 2.0f * lsb);
    }

    // Test edge cases
    {
        auto x = q16::from_float(0.0f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sqrt(0.0)", got, want, 2.0f * lsb);
    }

    // Test error handling for negative input
    {
        auto x = q16::from_float(-0.25f);
        auto result = x.sqrt();
        // Should return error value (most negative value for Q1.15 is -1.0)
        float got = result.to_float();
        float want = -1.0f;
        expect_near("sqrt(-0.25) error handling", got, want, 0.01f);
    }

    {
        auto x = q32::from_float(-1.0f);
        auto result = x.sqrt();
        // Should return error value (most negative value for Q3.29 is -4.0)
        float got = result.to_float();
        float want = -4.0f;
        expect_near("sqrt(-1.0) error handling Q32", got, want, 0.01f);
    }

    // Test Reference backend
    {
        auto x = q16_ref::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Reference sqrt(0.25)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16_ref::from_float(0.5f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Reference sqrt(0.5)", got, want, 2.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_sqrt_tests();
    return 0;
}
#endif
