#include "test_common.hpp"

namespace fp {
namespace test {

void run_antilogarithm_tests() {
    using q6_25 = q<6, 25, fp::test::Backend>;
    using q6_25_ref = q<6, 25, fp::test::Backend>;

    std::puts("\n--- Antilogarithm Tests ---");

    // Test antilog2 (2^x)
    // Input is interpreted as Q6.25, output is Q16.15
    {
        auto x = q6_25::from_float(1.0f);  // 2^1 = 2
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 2.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog2(1.0) = 2.0", got, want, 2.0f * lsb);
    }

    {
        auto x = q6_25::from_float(2.0f);  // 2^2 = 4
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 4.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog2(2.0) = 4.0", got, want, 2.0f * lsb);
    }

    {
        auto x = q6_25::from_float(-1.0f);  // 2^-1 = 0.5
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog2(-1.0) = 0.5", got, want, 2.0f * lsb);
    }

    // Test antilogn (e^x)
    {
        auto x = q6_25::from_float(1.0f);  // e^1 ≈ 2.718
        auto result = x.antilogn();
        float got = result.to_float();
        float want = std::exp(1.0f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilogn(1.0) ≈ 2.718", got, want, 2.0f * lsb);
    }

    {
        auto x = q6_25::from_float(0.0f);  // e^0 = 1
        auto result = x.antilogn();
        float got = result.to_float();
        float want = 1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilogn(0.0) = 1.0", got, want, 2.0f * lsb);
    }

    // Test antilog10 (10^x)
    {
        auto x = q6_25::from_float(1.0f);  // 10^1 = 10
        auto result = x.antilog10();
        float got = result.to_float();
        float want = 10.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog10(1.0) = 10.0", got, want, 2.0f * lsb);
    }

    {
        auto x = q6_25::from_float(0.0f);  // 10^0 = 1
        auto result = x.antilog10();
        float got = result.to_float();
        float want = 1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog10(0.0) = 1.0", got, want, 2.0f * lsb);
    }

    // Test with 32-bit input (priority 1 path for Xtensa)
    {
        auto x = q6_25::from_float(3.0f);  // 2^3 = 8
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog2(3.0) with Q32 (Xtensa prio-1)", got, want, 2.0f * lsb);
    }

    // Test Reference backend
    {
        auto x = q6_25_ref::from_float(2.0f);  // 2^2 = 4
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 4.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("antilog2(2.0) Reference backend", got, want, 2.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_antilogarithm_tests();
    return 0;
}
#endif
