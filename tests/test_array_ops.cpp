#include "test_common.hpp"

namespace fp {
namespace test {

void run_array_ops_tests() {
    using q16 = q<1, 15, Backend>;
    using q8  = q<1, 7, Backend>;
    using q16_ref = q<1, 15, ReferenceBackend>;

    std::puts("\n--- Array Min/Max Tests ---");

    // Test Option 1 API: FixedPointArray
    {
        int16_t data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(-0.25f).raw(),
            q16::from_float(0.75f).raw(),
            q16::from_float(0.1f).raw(),
            q16::from_float(-0.5f).raw()
        };

        fp::FixedPointArray<1, 15, Backend> arr(data, 5);

        auto min_val = arr.min();
        auto max_val = arr.max();

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.5f;
        float want_max = 0.75f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("FixedPointArray::min()", got_min, want_min, 2.0f * lsb);
        expect_near("FixedPointArray::max()", got_max, want_max, 2.0f * lsb);
    }

    // Test Option 3 API: Static methods on FixedPoint
    {
        int16_t data[] = {
            q16::from_float(0.3f).raw(),
            q16::from_float(0.9f).raw(),
            q16::from_float(-0.2f).raw(),
            q16::from_float(0.6f).raw()
        };

        auto min_val = q16::array_min(data, 4);
        auto max_val = q16::array_max(data, 4);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.2f;
        float want_max = 0.9f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("q16::array_min() static", got_min, want_min, 2.0f * lsb);
        expect_near("q16::array_max() static", got_max, want_max, 2.0f * lsb);
    }

    // Test with 8-bit arrays (priority 2 path for Xtensa)
    {
        int8_t data[] = {
            q8::from_float(0.5f).raw(),
            q8::from_float(0.25f).raw(),
            q8::from_float(-0.5f).raw(),
            q8::from_float(0.75f).raw()
        };

        auto min_val = q8::array_min(data, 4);
        auto max_val = q8::array_max(data, 4);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.5f;
        float want_max = 0.75f;

        const float lsb = 1.0f / static_cast<float>(1u << 7);
        expect_near("q8 array_min (Xtensa prio-2)", got_min, want_min, 2.0f * lsb);
        expect_near("q8 array_max (Xtensa prio-2)", got_max, want_max, 2.0f * lsb);
    }

    // Test with Reference backend
    {
        int16_t data[] = {
            q16_ref::from_float(0.2f).raw(),
            q16_ref::from_float(-0.7f).raw(),
            q16_ref::from_float(0.4f).raw()
        };

        auto min_val = q16_ref::array_min(data, 3);
        auto max_val = q16_ref::array_max(data, 3);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.7f;
        float want_max = 0.4f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Reference backend array_min", got_min, want_min, 2.0f * lsb);
        expect_near("Reference backend array_max", got_max, want_max, 2.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_array_ops_tests();
    return 0;
}
#endif
