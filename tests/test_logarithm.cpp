#include "test_common.hpp"

namespace fp {
namespace test {

void run_logarithm_tests() {
    using q16 = q<1, 15, fp::test::Backend>;
    using q32 = q<3, 29, fp::test::Backend>;

    std::puts("\n--- Logarithm Tests ---");

    // Basic log2 tests
    {
        auto x = q16::from_float(0.5f);  // log2(0.5) = -1
        auto result = x.log2();
        float got = result.to_float();
        float want = -1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        expect_near("log2(0.5)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16::from_float(0.25f);  // log2(0.25) = -2
        auto result = x.log2();
        float got = result.to_float();
        float want = -2.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        expect_near("log2(0.25)", got, want, 2.0f * lsb);
    }

    // Basic logn test
    {
        auto x = q16::from_float(0.5f);  // logn(0.5) ≈ -0.693147
        auto result = x.logn();
        float got = result.to_float();
        float want = std::log(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        expect_near("logn(0.5)", got, want, 2.0f * lsb);
    }

    // Basic log10 test
    {
        auto x = q16::from_float(0.1f);  // log10(0.1) = -1
        auto result = x.log10();
        float got = result.to_float();
        float want = -1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        expect_near("log10(0.1)", got, want, 1000.0f * lsb);  // More lenient epsilon
    }

    // Test with 32-bit input
    {
        auto x = q32::from_float(0.5f);
        auto result = x.log2();
        float got = result.to_float();
        float want = -1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        expect_near("log2(0.5) with Q32", got, want, 2.0f * lsb);
    }

    // Test error handling for zero/negative
    {
        auto x = q16::from_float(-0.5f);  // negative input
        auto result = x.log2();
        // Result should be most negative value (0x80000000)
        // which as Q6.25 is -64.0
        float got = result.to_float();
        float want = -64.0f;  // Approximate expected error value
        expect_near("log2(-0.5) error handling", got, want, 1.0f);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_logarithm_tests();
    return 0;
}
#endif
