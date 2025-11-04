/**
 * Fixed-Point Library Educational Tutorial
 *
 * This example demonstrates:
 * 1. Basic fixed-point operations
 * 2. Transparent mixed-precision handling
 * 3. Scalar vs. vector (array) operations
 * 4. How the library handles requantization automatically
 */

#include "../fp.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

void print_header(const char* title) {
    std::cout << "\n========================================\n";
    std::cout << title << "\n";
    std::cout << "========================================\n";
}

void print_value(const char* label, float value, const char* unit = "") {
    std::cout << std::setw(30) << std::left << label
              << ": " << std::fixed << std::setprecision(6)
              << value << " " << unit << "\n";
}

int main() {
    using namespace fp;

    print_header("PART 1: Basic Operations - Same Precision");
    {
        // Q1.15 format: 1 integer bit, 15 fractional bits
        // Range: [-1.0, 1.0), Resolution: 1/32768 ≈ 0.000031
        using Q1_15 = q<1, 15>;

        auto a = Q1_15::from_float(0.5f);
        auto b = Q1_15::from_float(0.75f);

        // Multiplication (result defaults to same format as left operand)
        auto product = a * b;
        print_value("0.5 × 0.75", product.to_float());
        print_value("Expected", 0.5f * 0.75f);
        print_value("Error", std::abs(product.to_float() - 0.375f));

        // Division (result saturates in Q1.15 which can only hold [-1, 1))
        auto quotient = b / a;  // Result would be 1.5, but saturates to ~1.0
        print_value("0.75 ÷ 0.5 (saturated in Q1.15)", quotient.to_float());

        // Use larger output format to avoid saturation
        using Q3_13 = q<3, 13>;  // Can hold [-4, 4)
        auto quotient_unsaturated = div_as<3, 13>(b, a);
        print_value("0.75 ÷ 0.5 (using Q3.13)", quotient_unsaturated.to_float());
        print_value("Expected", 0.75f / 0.5f);
    }

    print_header("PART 2: Mixed Precision - Transparent Requantization");
    {
        // Different Q formats can be mixed seamlessly!
        // The library handles all requantization and shifting automatically.

        using Q4_12 = q<4, 12>;  // Range: [-8, 8), Resolution: 1/4096
        using Q1_15 = q<1, 15>;  // Range: [-1, 1), Resolution: 1/32768
        using Q3_13 = q<3, 13>;  // Range: [-4, 4), Resolution: 1/8192

        auto large_val = Q4_12::from_float(3.5f);   // Needs 4 integer bits
        auto small_val = Q1_15::from_float(0.25f);  // Fits in 1 integer bit

        std::cout << "Multiplying Q4.12 (3.5) × Q1.15 (0.25):\n";

        // Result defaults to left operand's format (Q4.12)
        auto result1 = large_val * small_val;
        print_value("  Result (Q4.12 format)", result1.to_float());

        // We can explicitly specify output format using mul_as<>
        auto result2 = mul_as<3, 13>(large_val, small_val);
        print_value("  Result (Q3.13 format)", result2.to_float());

        // Both give same answer, just different internal representation
        std::cout << "\nKey Point: No manual shifting or requantization needed!\n";
        std::cout << "The library handles all the fixed-point arithmetic.\n";
    }

    print_header("PART 3: Scalar Operations (Element-wise)");
    {
        using Q1_15 = q<1, 15>;

        // Scalar operations work on individual values
        auto x = Q1_15::from_float(0.5f);

        std::cout << "Starting with x = 0.5:\n";

        // Logarithm (returns Q6.25 format for extended range)
        auto log_x = x.log2();
        print_value("  log2(x)", log_x.to_float());
        print_value("  Expected", std::log2(0.5f));

        // Square root (returns same format as input)
        auto sqrt_x = x.sqrt();
        print_value("  sqrt(x)", sqrt_x.to_float());
        print_value("  Expected", std::sqrt(0.5f));

        // Power operation
        auto exp = Q1_15::from_float(0.5f);
        auto pow_x = x.pow(exp);
        print_value("  x^0.5 (same as sqrt)", pow_x.to_float());

        // Trigonometric functions (input in radians)
        auto angle = Q1_15::from_float(0.7854f);  // π/4
        auto sin_result = angle.sin();
        print_value("  sin(π/4)", sin_result.to_float());
        print_value("  Expected", std::sin(0.7854f));

        // Activation functions
        auto sigmoid_result = x.sigmoid();
        print_value("  sigmoid(x)", sigmoid_result.to_float());
    }

    print_header("PART 4: Vector Operations (Array-based)");
    {
        using Q1_15 = q<1, 15>;

        // Create array of fixed-point values (stored as int16_t internally)
        int16_t data[] = {
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(-0.25f).raw(),
            Q1_15::from_float(0.75f).raw(),
            Q1_15::from_float(0.1f).raw(),
            Q1_15::from_float(-0.8f).raw()
        };

        std::cout << "Array: [0.5, -0.25, 0.75, 0.1, -0.8]\n\n";

        // Option 1: Using FixedPointArray wrapper
        {
            FixedPointArray<1, 15> arr(data, 5);

            auto min_val = arr.min();
            auto max_val = arr.max();
            auto sum_val = arr.sum();

            print_value("  Min", min_val.to_float());
            print_value("  Max", max_val.to_float());
            print_value("  Sum", sum_val.to_float());
        }

        // Option 2: Using static methods
        {
            auto min_val = Q1_15::array_min(data, 5);
            auto max_val = Q1_15::array_max(data, 5);

            std::cout << "\nUsing static methods (alternative API):\n";
            print_value("  Min", min_val.to_float());
            print_value("  Max", max_val.to_float());
        }

        // In-place array operations
        {
            int16_t arr_copy[5];
            std::copy(data, data + 5, arr_copy);

            FixedPointArray<1, 15> arr(arr_copy, 5);

            // Scale all elements by 2.0
            auto scale_factor = Q1_15::from_float(0.5f);  // Will shift left by 1 bit (×2)
            arr.scale(scale_factor);

            std::cout << "\nAfter scaling by 2.0:\n";
            std::cout << "  Scaled array: [";
            for (int i = 0; i < 5; i++) {
                std::cout << Q1_15(arr_copy[i]).to_float();
                if (i < 4) std::cout << ", ";
            }
            std::cout << "]\n";
        }

        // Dot product (vector operation)
        {
            int16_t vec1[] = {
                Q1_15::from_float(0.2f).raw(),
                Q1_15::from_float(0.3f).raw(),
                Q1_15::from_float(0.1f).raw()
            };

            int16_t vec2[] = {
                Q1_15::from_float(0.4f).raw(),
                Q1_15::from_float(-0.2f).raw(),
                Q1_15::from_float(0.5f).raw()
            };

            FixedPointArray<1, 15> array1(vec1, 3);
            FixedPointArray<1, 15> array2(vec2, 3);

            auto dot = array1.dot_product(array2);

            std::cout << "\nDot product (NOTE: raw integer multiply, may have precision loss):\n";
            std::cout << "  [0.2, 0.3, 0.1] · [0.4, -0.2, 0.5]\n";
            print_value("  Result", dot.to_float());
            print_value("  Expected", 0.2f*0.4f + 0.3f*(-0.2f) + 0.1f*0.5f);
        }
    }

    print_header("PART 5: Complex Example - Neural Network Layer");
    {
        using Q1_15 = q<1, 15>;

        std::cout << "Simple neural network: weighted sum + activation\n\n";

        // Input features [0.3, 0.5, 0.2]
        int16_t inputs[] = {
            Q1_15::from_float(0.3f).raw(),
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(0.2f).raw()
        };

        // Weights [0.4, -0.3, 0.6]
        int16_t weights[] = {
            Q1_15::from_float(0.4f).raw(),
            Q1_15::from_float(-0.3f).raw(),
            Q1_15::from_float(0.6f).raw()
        };

        FixedPointArray<1, 15> input_array(inputs, 3);
        FixedPointArray<1, 15> weight_array(weights, 3);

        // Step 1: Compute weighted sum (dot product)
        auto weighted_sum = input_array.dot_product(weight_array);
        print_value("  Weighted sum", weighted_sum.to_float());
        print_value("  Expected", 0.3f*0.4f + 0.5f*(-0.3f) + 0.2f*0.6f);

        // Step 2: Apply sigmoid activation
        auto output = weighted_sum.sigmoid();
        print_value("  After sigmoid", output.to_float());

        std::cout << "\nKey Point: Complex operations composed from primitives!\n";
    }

    print_header("PART 6: Backend Selection (Reference vs Xtensa)");
    {
        // The library supports multiple backends:
        // - ReferenceBackend: Portable C++ implementation
        // - XtensaBackend: Optimized for Cadence Xtensa HiFi DSP

        using Q16_ref = q<1, 15, ReferenceBackend>;

        auto a_ref = Q16_ref::from_float(0.5f);
        auto b_ref = Q16_ref::from_float(0.75f);
        auto result_ref = a_ref * b_ref;

#ifdef __XTENSA__
        using Q16_xtensa = q<1, 15, XtensaBackend>;
        auto a_xtensa = Q16_xtensa::from_float(0.5f);
        auto b_xtensa = Q16_xtensa::from_float(0.75f);
        auto result_xtensa = a_xtensa * b_xtensa;

        std::cout << "Both backends produce same results:\n";
        print_value("  Reference backend", result_ref.to_float());
        print_value("  Xtensa backend", result_xtensa.to_float());
#else
        std::cout << "Reference backend result:\n";
        print_value("  Value", result_ref.to_float());
#endif

        std::cout << "\nKey Point: API is identical, backend selection\n";
        std::cout << "only affects performance, not functionality!\n";
    }

    print_header("SUMMARY");
    std::cout << R"(
This library provides:

1. ✓ Type-safe fixed-point arithmetic
2. ✓ Automatic mixed-precision handling
3. ✓ No manual shifting/requantization needed
4. ✓ Both scalar and vector operations
5. ✓ Multiple backend support (portable & optimized)
6. ✓ Clean, intuitive API

Common patterns:
- Scalar: auto result = a.operation(b);
- Vector: FixedPointArray<I,F> arr(data, len);
- Mixed precision: auto result = mul_as<OUT_I, OUT_F>(a, b);

For more details, see fp.hpp and test files in tests/
)";

    return 0;
}
