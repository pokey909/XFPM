#include "test_common.hpp"

namespace fp {
namespace test {

void run_vector_ops_tests() {
    using q16 = q<1, 15>;  // Q1.15 format (16-bit total)
    using q8  = q<1, 7>;   // Q1.7 format (8-bit total)

    std::puts("\n--- Vector Element-wise Tests ---");

    // Test element-wise multiply (Q1.15 format: values must be in [-1.0, ~1.0))
    {
        int16_t arr1_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.25f).raw(),
            q16::from_float(0.75f).raw(),
            q16::from_float(-0.5f).raw()
        };
        int16_t arr2_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(-0.5f).raw()
        };
        int16_t output_data[4];

        fp::q_array<1, 15> arr1(arr1_data, 4);
        fp::q_array<1, 15> arr2(arr2_data, 4);
        fp::q_array<1, 15> output(output_data, 4);

        arr1.elemult(arr2, output);

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("elemult [0]: 0.5*0.5", output[0].to_float(), 0.25f, 2.0f * lsb);
        expect_near("elemult [1]: 0.25*0.5", output[1].to_float(), 0.125f, 2.0f * lsb);
        expect_near("elemult [2]: 0.75*0.5", output[2].to_float(), 0.375f, 2.0f * lsb);
        expect_near("elemult [3]: -0.5*-0.5", output[3].to_float(), 0.25f, 2.0f * lsb);
    }

    // Test element-wise add
    {
        int16_t arr1_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.25f).raw(),
            q16::from_float(-0.5f).raw()
        };
        int16_t arr2_data[] = {
            q16::from_float(0.25f).raw(),
            q16::from_float(0.75f).raw(),
            q16::from_float(0.5f).raw()
        };
        int16_t output_data[3];

        fp::q_array<1, 15> arr1(arr1_data, 3);
        fp::q_array<1, 15> arr2(arr2_data, 3);
        fp::q_array<1, 15> output(output_data, 3);

        arr1.add(arr2, output);

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("add [0]: 0.5+0.25", output[0].to_float(), 0.75f, 2.0f * lsb);
        expect_near("add [1]: 0.25+0.75", output[1].to_float(), 1.0f, 2.0f * lsb);
        expect_near("add [2]: -0.5+0.5", output[2].to_float(), 0.0f, 2.0f * lsb);
    }

    // Test element-wise subtract
    {
        int16_t arr1_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(1.0f).raw(),
            q16::from_float(0.25f).raw()
        };
        int16_t arr2_data[] = {
            q16::from_float(0.25f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(0.75f).raw()
        };
        int16_t output_data[3];

        fp::q_array<1, 15> arr1(arr1_data, 3);
        fp::q_array<1, 15> arr2(arr2_data, 3);
        fp::q_array<1, 15> output(output_data, 3);

        arr1.sub(arr2, output);

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sub [0]: 0.5-0.25", output[0].to_float(), 0.25f, 2.0f * lsb);
        expect_near("sub [1]: 1.0-0.5", output[1].to_float(), 0.5f, 2.0f * lsb);
        expect_near("sub [2]: 0.25-0.75", output[2].to_float(), -0.5f, 2.0f * lsb);
    }

    std::puts("\n--- Vector Statistical Tests ---");

    // Test mean
    {
        int16_t arr_data[] = {
            q16::from_float(0.0f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(1.0f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 3);
        auto mean = arr.mean();

        // Mean should be (0 + 0.5 + 1.0) / 3 = 0.5
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("mean of [0, 0.5, 1.0]", mean.to_float(), 0.5f, 2.0f * lsb);
    }

    // Test mean with negative values
    {
        int16_t arr_data[] = {
            q16::from_float(-1.0f).raw(),
            q16::from_float(-0.5f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(1.0f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 4);
        auto mean = arr.mean();

        // Mean should be (-1.0 + -0.5 + 0.5 + 1.0) / 4 = 0.0
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("mean of [-1, -0.5, 0.5, 1]", mean.to_float(), 0.0f, 2.0f * lsb);
    }

    // Test RMS
    {
        int16_t arr_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.5f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 2);
        auto rms = arr.rms();

        // RMS = sqrt((0.5^2 + 0.5^2) / 2) = sqrt(0.5/2) = sqrt(0.25) = 0.5
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("rms of [0.5, 0.5]", rms.to_float(), 0.5f, 4.0f * lsb);
    }

    // Test variance
    {
        int16_t arr_data[] = {
            q16::from_float(0.0f).raw(),
            q16::from_float(0.5f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 2);
        auto variance = arr.variance();

        // Mean = 0.25
        // Variance = ((0-0.25)^2 + (0.5-0.25)^2) / 2 = (0.0625 + 0.0625) / 2 = 0.0625
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("variance of [0, 0.5]", variance.to_float(), 0.0625f, 4.0f * lsb);
    }

    // Test standard deviation
    {
        int16_t arr_data[] = {
            q16::from_float(0.0f).raw(),
            q16::from_float(0.5f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 2);
        auto stddev = arr.stddev();

        // Stddev = sqrt(variance) = sqrt(0.0625) = 0.25
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("stddev of [0, 0.5]", stddev.to_float(), 0.25f, 4.0f * lsb);
    }

    // Test with Q8 format (8-bit path)
    {
        int8_t arr_data[] = {
            q8::from_float(0.5f).raw(),
            q8::from_float(0.5f).raw()
        };

        fp::q_array<1, 7> arr(arr_data, 2);
        auto mean = arr.mean();

        // Mean should be (0.5 + 0.5) / 2 = 0.5
        const float lsb = 1.0f / static_cast<float>(1u << 7);
        expect_near("Q8 mean", mean.to_float(), 0.5f, 2.0f * lsb);
    }

    // Test Backend alias explicitly (length must be multiple of 4)
    {
        using q16_backend = fp::q<1, 15, test::Backend>;
        int16_t arr1_data[] = {
            q16_backend::from_float(0.5f).raw(),
            q16_backend::from_float(0.75f).raw(),
            q16_backend::from_float(0.25f).raw(),
            q16_backend::from_float(0.5f).raw()
        };
        int16_t arr2_data[] = {
            q16_backend::from_float(0.5f).raw(),
            q16_backend::from_float(0.5f).raw(),
            q16_backend::from_float(0.5f).raw(),
            q16_backend::from_float(0.25f).raw()
        };
        int16_t output_data[4] = {0, 0, 0, 0};

        fp::q_array<1, 15, test::Backend> arr1(arr1_data, 4);
        fp::q_array<1, 15, test::Backend> arr2(arr2_data, 4);
        fp::q_array<1, 15, test::Backend> output(output_data, 4);

        arr1.elemult(arr2, output);

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Backend elemult [0]", output[0].to_float(), 0.25f, 2.0f * lsb);
        expect_near("Backend elemult [1]", output[1].to_float(), 0.375f, 2.0f * lsb);
        expect_near("Backend elemult [2]", output[2].to_float(), 0.125f, 2.0f * lsb);
        expect_near("Backend elemult [3]", output[3].to_float(), 0.125f, 2.0f * lsb);
    }

    // Test Reference backend explicitly
    {
        using q16_ref = fp::q<1, 15, fp::ReferenceBackend>;
        int16_t arr_data[] = {
            q16_ref::from_float(0.5f).raw(),
            q16_ref::from_float(0.75f).raw()
        };

        fp::q_array<1, 15, fp::ReferenceBackend> arr(arr_data, 2);
        auto mean = arr.mean();

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("Reference mean", mean.to_float(), 0.625f, 2.0f * lsb);
    }

    // Test sum operation
    {
        int16_t arr_data[] = {
            q16::from_float(0.25f).raw(),
            q16::from_float(0.5f).raw(),
            q16::from_float(0.125f).raw()
        };

        fp::q_array<1, 15> arr(arr_data, 3);
        auto sum = arr.sum();

        // Sum should be 0.25 + 0.5 + 0.125 = 0.875
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("sum of [0.25, 0.5, 0.125]", sum.to_float(), 0.875f, 2.0f * lsb);
    }

    // Test dot product
    {
        int16_t arr1_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.25f).raw()
        };
        int16_t arr2_data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(0.5f).raw()
        };

        fp::q_array<1, 15> arr1(arr1_data, 2);
        fp::q_array<1, 15> arr2(arr2_data, 2);
        auto dot = arr1.dot_product(arr2);

        // Dot product = 0.5*0.5 + 0.25*0.5 = 0.25 + 0.125 = 0.375
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        expect_near("dot product [0.5,0.25]·[0.5,0.5]", dot.to_float(), 0.375f, 4.0f * lsb);
    }
}

} // namespace test
} // namespace fp

// Main function for standalone execution (not used when building full test suite)
#ifndef FP_TEST_SUITE
int main() {
    fp::test::run_vector_ops_tests();
    return 0;
}
#endif
