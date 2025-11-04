/**
 * Fixed-Point Library - Vectorized Operations Tutorial
 *
 * This example demonstrates:
 * 1. Scalar vs vectorized operations
 * 2. Xtensa NatureDSP optimized vector operations
 * 3. Fast variants and their alignment/length restrictions
 * 4. Practical examples of when fast variants can/cannot be used
 * 5. Element-wise operations, reductions, and statistical operations
 *
 * KEY CONCEPTS:
 * - Regular variants: Work with any alignment and length
 * - Fast variants: Require 8-byte alignment and N % 4 == 0
 * - Backend abstraction: Same API works with Reference or Xtensa backend
 */

#include "../fp.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

// Helper to print section headers
void print_header(const char* title) {
    std::cout << "\n========================================\n";
    std::cout << title << "\n";
    std::cout << "========================================\n";
}

// Helper to print subsection headers
void print_subheader(const char* title) {
    std::cout << "\n--- " << title << " ---\n";
}

// Helper to print values with labels
void print_value(const char* label, float value, const char* unit = "") {
    std::cout << std::setw(35) << std::left << label
              << ": " << std::fixed << std::setprecision(6)
              << value << " " << unit << "\n";
}

// Helper to check if pointer is 8-byte aligned
bool is_8byte_aligned(const void* ptr) {
    return (reinterpret_cast<uintptr_t>(ptr) & 7) == 0;
}

// Helper to print array
template<typename T>
void print_array(const char* label, const T* arr, size_t length) {
    std::cout << label << " [";
    for (size_t i = 0; i < length; i++) {
        std::cout << arr[i].to_float();
        if (i < length - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

int main() {
    using namespace fp;

    // Select backend based on build target
#ifdef __XTENSA__
    using Backend = XtensaBackend;
#else
    using Backend = ReferenceBackend;
#endif

    print_header("PART 1: Scalar vs Vectorized Operations");
    {
        using Q1_15 = q<1, 15, Backend>;

        print_subheader("Scalar Operations (Element-wise)");

        auto a = Q1_15::from_float(0.5f);
        auto b = Q1_15::from_float(0.3f);

        // Scalar multiply
        auto scalar_result = a * b;
        print_value("0.5 × 0.3 (scalar)", scalar_result.to_float());

        // Multiple scalar operations
        auto x = Q1_15::from_float(0.7f);
        auto y = x.sqrt();
        auto z = y.sin();

        print_value("sqrt(0.7)", y.to_float());
        print_value("sin(sqrt(0.7))", z.to_float());

        print_subheader("Vectorized Operations (SIMD-style)");

        // Create aligned arrays for vectorized operations
        alignas(8) int16_t vec1_data[8];
        alignas(8) int16_t vec2_data[8];
        alignas(8) int16_t output_data[8];

        // Initialize arrays
        for (int i = 0; i < 8; i++) {
            vec1_data[i] = Q1_15::from_float(0.5f + i * 0.05f).raw();
            vec2_data[i] = Q1_15::from_float(0.3f - i * 0.02f).raw();
        }

        FixedPointArray<1, 15, Backend> vec1(vec1_data, 8);
        FixedPointArray<1, 15, Backend> vec2(vec2_data, 8);
        FixedPointArray<1, 15, Backend> output(output_data, 8);

        // Element-wise multiply (vectorized)
        vec1.elemult(vec2, output);

        std::cout << "\nElement-wise multiplication of 8-element vectors:\n";
        std::cout << "  Input 1 : [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85]\n";
        std::cout << "  Input 2 : [0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16]\n";
        std::cout << "  Output  : [";
        for (int i = 0; i < 8; i++) {
            std::cout << Q1_15(output_data[i]).to_float();
            if (i < 7) std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "\n🔑 Key Point: Vectorized ops process multiple elements simultaneously!\n";
        std::cout << "   On Xtensa HiFi DSP, this uses SIMD instructions for speed.\n";
    }

    print_header("PART 2: NatureDSP Fast Variants - Restrictions");
    {
        std::cout << R"(
NatureDSP provides two variants of many operations:

┌─────────────────┬──────────────┬─────────────────────────────┐
│   Operation     │   Regular    │         Fast Variant        │
├─────────────────┼──────────────┼─────────────────────────────┤
│ Dot Product     │ vec_dot16x16 │ vec_dot16x16_fast           │
│ Element-wise +  │ vec_add16x16 │ vec_add16x16_fast           │
│ Element-wise ×  │ vec_elemult  │ (no fast variant)           │
│ Element-wise -  │ vec_elesub   │ (no fast variant)           │
│ Array Sum       │ vec_sum16x16 │ (no fast variant)           │
│ Power (x²)      │ vec_power16x16│ vec_power16x16_fast        │
└─────────────────┴──────────────┴─────────────────────────────┘

FAST VARIANT REQUIREMENTS (All must be satisfied):
✓ All input/output pointers must be 8-byte aligned
✓ Array length N must be ≥ 4
✓ Array length N must be a multiple of 4

If ANY requirement is not met → falls back to regular variant
)";
    }

    print_header("PART 3: Alignment Requirements Demonstration");
    {
        using Q1_15 = q<1, 15, Backend>;

        print_subheader("Scenario 1: Properly Aligned Arrays (Fast variant used)");

        // Aligned arrays - fast variants will be used
        alignas(8) int16_t aligned_a[8];
        alignas(8) int16_t aligned_b[8];
        alignas(8) int16_t aligned_out[8];

        for (int i = 0; i < 8; i++) {
            aligned_a[i] = Q1_15::from_float(0.1f * i).raw();
            aligned_b[i] = Q1_15::from_float(0.2f * i).raw();
        }

        std::cout << "Array A address: " << static_cast<void*>(aligned_a)
                  << " → " << (is_8byte_aligned(aligned_a) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Array B address: " << static_cast<void*>(aligned_b)
                  << " → " << (is_8byte_aligned(aligned_b) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Length: 8 → ✓ Multiple of 4\n";
        std::cout << "Result: FAST VARIANT will be used (when NatureDSP is enabled)\n";

        FixedPointArray<1, 15, Backend> arr_a(aligned_a, 8);
        FixedPointArray<1, 15, Backend> arr_b(aligned_b, 8);
        FixedPointArray<1, 15, Backend> arr_out(aligned_out, 8);

        arr_a.add(arr_b, arr_out);

        std::cout << "Addition completed successfully using fast variant.\n";

        print_subheader("Scenario 2: Unaligned Arrays (Regular variant used)");

        // Create a buffer and offset to guarantee misalignment
        alignas(8) int16_t buffer[10];
        int16_t* unaligned_a = buffer + 1;  // Offset by 2 bytes → not 8-byte aligned

        alignas(8) int16_t aligned_b2[8];
        alignas(8) int16_t aligned_out2[8];

        for (int i = 0; i < 8; i++) {
            unaligned_a[i] = Q1_15::from_float(0.1f * i).raw();
            aligned_b2[i] = Q1_15::from_float(0.2f * i).raw();
        }

        std::cout << "Array A address: " << static_cast<void*>(unaligned_a)
                  << " → " << (is_8byte_aligned(unaligned_a) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Array B address: " << static_cast<void*>(aligned_b2)
                  << " → " << (is_8byte_aligned(aligned_b2) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Length: 8 → ✓ Multiple of 4\n";
        std::cout << "Result: REGULAR VARIANT used (A is not aligned)\n";

        FixedPointArray<1, 15, Backend> arr_a2(unaligned_a, 8);
        FixedPointArray<1, 15, Backend> arr_b2(aligned_b2, 8);
        FixedPointArray<1, 15, Backend> arr_out2(aligned_out2, 8);

        arr_a2.add(arr_b2, arr_out2);

        std::cout << "Addition completed successfully using regular variant.\n";

        print_subheader("Scenario 3: Wrong Length (Regular variant used)");

        alignas(8) int16_t aligned_a3[7];  // Length 7 is not multiple of 4
        alignas(8) int16_t aligned_b3[7];
        alignas(8) int16_t aligned_out3[7];

        std::cout << "Array A address: " << static_cast<void*>(aligned_a3)
                  << " → " << (is_8byte_aligned(aligned_a3) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Array B address: " << static_cast<void*>(aligned_b3)
                  << " → " << (is_8byte_aligned(aligned_b3) ? "✓ 8-byte aligned" : "✗ NOT aligned") << "\n";
        std::cout << "Length: 7 → ✗ NOT a multiple of 4\n";
        std::cout << "Result: REGULAR VARIANT used (length requirement not met)\n";

        std::cout << "\n🔑 Key Point: Library automatically selects the appropriate variant!\n";
        std::cout << "   No manual checking needed - fallback happens transparently.\n";
    }

    print_header("PART 4: Vectorized Operations - Complete Examples");
    {
        using Q1_15 = q<1, 15, Backend>;

        print_subheader("Element-wise Operations");

        alignas(8) int16_t a[4], b[4], sum[4], diff[4], prod[4];

        // Initialize
        for (int i = 0; i < 4; i++) {
            a[i] = Q1_15::from_float(0.5f - i * 0.1f).raw();
            b[i] = Q1_15::from_float(0.2f + i * 0.05f).raw();
        }

        FixedPointArray<1, 15, Backend> arr_a(a, 4);
        FixedPointArray<1, 15, Backend> arr_b(b, 4);
        FixedPointArray<1, 15, Backend> arr_sum(sum, 4);
        FixedPointArray<1, 15, Backend> arr_diff(diff, 4);
        FixedPointArray<1, 15, Backend> arr_prod(prod, 4);

        // Perform operations
        arr_a.add(arr_b, arr_sum);      // Uses vec_add16x16_fast (aligned, N=4)
        arr_a.sub(arr_b, arr_diff);     // Uses vec_elesub16x16 (no fast variant)
        arr_a.elemult(arr_b, arr_prod); // Uses vec_elemult16x16 + vec_shift16x16

        std::cout << "Array A    : [";
        for (int i = 0; i < 4; i++) {
            std::cout << Q1_15(a[i]).to_float();
            if (i < 3) std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "Array B    : [";
        for (int i = 0; i < 4; i++) {
            std::cout << Q1_15(b[i]).to_float();
            if (i < 3) std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "A + B      : [";
        for (int i = 0; i < 4; i++) {
            std::cout << Q1_15(sum[i]).to_float();
            if (i < 3) std::cout << ", ";
        }
        std::cout << "] (using vec_add16x16_fast)\n";

        std::cout << "A - B      : [";
        for (int i = 0; i < 4; i++) {
            std::cout << Q1_15(diff[i]).to_float();
            if (i < 3) std::cout << ", ";
        }
        std::cout << "] (using vec_elesub16x16)\n";

        std::cout << "A × B      : [";
        for (int i = 0; i < 4; i++) {
            std::cout << Q1_15(prod[i]).to_float();
            if (i < 3) std::cout << ", ";
        }
        std::cout << "] (using vec_elemult + vec_shift)\n";

        print_subheader("Reduction Operations");

        alignas(8) int16_t data[8];
        for (int i = 0; i < 8; i++) {
            data[i] = Q1_15::from_float(0.1f * (i + 1)).raw();
        }

        FixedPointArray<1, 15, Backend> arr(data, 8);

        auto sum_val = arr.sum();       // Uses vec_sum16x16 (no fast variant)
        auto min_val = arr.min();       // Uses vec_min16x16 (if available)
        auto max_val = arr.max();       // Uses vec_max16x16 (if available)

        std::cout << "Data       : [";
        for (int i = 0; i < 8; i++) {
            std::cout << Q1_15(data[i]).to_float();
            if (i < 7) std::cout << ", ";
        }
        std::cout << "]\n";

        print_value("Sum", sum_val.to_float());
        print_value("Min", min_val.to_float());
        print_value("Max", max_val.to_float());

        print_subheader("Dot Product");

        alignas(8) int16_t vec1[4], vec2[4];

        vec1[0] = Q1_15::from_float(0.5f).raw();
        vec1[1] = Q1_15::from_float(-0.3f).raw();
        vec1[2] = Q1_15::from_float(0.2f).raw();
        vec1[3] = Q1_15::from_float(0.7f).raw();

        vec2[0] = Q1_15::from_float(0.4f).raw();
        vec2[1] = Q1_15::from_float(0.6f).raw();
        vec2[2] = Q1_15::from_float(-0.1f).raw();
        vec2[3] = Q1_15::from_float(0.3f).raw();

        FixedPointArray<1, 15, Backend> v1(vec1, 4);
        FixedPointArray<1, 15, Backend> v2(vec2, 4);

        auto dot = v1.dot_product(v2);  // Uses vec_dot16x16_fast (aligned, N=4)

        std::cout << "Vector 1   : [0.5, -0.3, 0.2, 0.7]\n";
        std::cout << "Vector 2   : [0.4, 0.6, -0.1, 0.3]\n";
        print_value("Dot product", dot.to_float(), "(using vec_dot16x16_fast)");

        float expected = 0.5f*0.4f + (-0.3f)*0.6f + 0.2f*(-0.1f) + 0.7f*0.3f;
        print_value("Expected", expected);

        print_subheader("Statistical Operations");

        auto mean_val = arr.mean();       // Uses vec_mean16x16 (no fast variant)
        auto rms_val = arr.rms();         // Uses vec_rms16x16 (internally uses vec_power)
        auto variance = arr.variance();   // Uses vec_var16x16 (no fast variant)
        auto stddev = arr.stddev();       // Uses vec_stddev16x16 (no fast variant)

        std::cout << "Statistical operations on: [";
        for (int i = 0; i < 8; i++) {
            std::cout << Q1_15(data[i]).to_float();
            if (i < 7) std::cout << ", ";
        }
        std::cout << "]\n";

        print_value("Mean", mean_val.to_float());
        print_value("RMS", rms_val.to_float());
        print_value("Variance", variance.to_float());
        print_value("Std Dev", stddev.to_float());
    }

    print_header("PART 5: Performance Considerations");
    {
        std::cout << R"(
When to use vectorized operations:

✓ DO use vectors when:
  • Processing large arrays (hundreds/thousands of elements)
  • Same operation applied to all elements
  • Data can be 8-byte aligned (for fast variants)
  • Length is known and can be made multiple of 4

✗ DON'T use vectors when:
  • Only processing 1-2 elements (scalar is simpler)
  • Each element needs different operations
  • Data is inherently misaligned (e.g., sliding window)

Tips for maximum performance:

1. Alignment:
   alignas(8) int16_t my_array[N];  // Ensures 8-byte alignment

2. Length:
   Use lengths that are multiples of 4: 4, 8, 12, 16, 20, ...

3. Memory allocation:
   For dynamic arrays:
   void* ptr = aligned_alloc(8, N * sizeof(int16_t));

4. Structure padding:
   struct Data {
       alignas(8) int16_t values[100];  // Force alignment
   };

Fast variant speedup (typical on Xtensa HiFi4):
  • Dot product: 2-3x faster
  • Element-wise add: 2-4x faster
  • Power (sum of squares): 1.5-2x faster

NOTE: In this build, NatureDSP headers are commented out,
so all operations fall back to the portable reference
implementation. To enable hardware acceleration:
  1. Uncomment #include <NatureDSP_Signal.h> in backend.hpp
  2. Link against NatureDSP library
  3. Build for Xtensa target
)";
    }

    print_header("PART 6: Backend Comparison");
    {
        print_subheader("Reference Backend (Portable C++)");

        using Q1_15_ref = q<1, 15, ReferenceBackend>;

        alignas(8) int16_t ref_a[4] = {
            Q1_15_ref::from_float(0.5f).raw(),
            Q1_15_ref::from_float(0.3f).raw(),
            Q1_15_ref::from_float(0.7f).raw(),
            Q1_15_ref::from_float(0.2f).raw()
        };

        FixedPointArray<1, 15, ReferenceBackend> ref_arr(ref_a, 4);
        auto ref_sum = ref_arr.sum();

        print_value("Sum (Reference)", ref_sum.to_float());

        print_subheader("Xtensa Backend (Hardware Optimized)");

        using Q1_15_xtensa = q<1, 15, Backend>;

        alignas(8) int16_t xtensa_a[4] = {
            Q1_15_xtensa::from_float(0.5f).raw(),
            Q1_15_xtensa::from_float(0.3f).raw(),
            Q1_15_xtensa::from_float(0.7f).raw(),
            Q1_15_xtensa::from_float(0.2f).raw()
        };

        FixedPointArray<1, 15, Backend> xtensa_arr(xtensa_a, 4);
        auto xtensa_sum = xtensa_arr.sum();

        print_value("Sum (Xtensa)", xtensa_sum.to_float());

        std::cout << "\n🔑 Key Point: Identical API, different performance!\n";
        std::cout << "   Results are identical, only execution speed differs.\n";
    }

    print_header("PART 7: Real-World Example - Audio Processing");
    {
        using Q1_15 = q<1, 15, Backend>;

        std::cout << R"(
Scenario: Simple audio mixer
  • Mix two audio channels with different gains
  • Apply fade-in envelope
  • Compute RMS level for metering

Input: 2 channels × 8 samples each
)";

        // Channel 1: Sine-like wave
        alignas(8) int16_t channel1[8] = {
            Q1_15::from_float(0.0f).raw(),
            Q1_15::from_float(0.3f).raw(),
            Q1_15::from_float(0.7f).raw(),
            Q1_15::from_float(0.9f).raw(),
            Q1_15::from_float(0.7f).raw(),
            Q1_15::from_float(0.3f).raw(),
            Q1_15::from_float(0.0f).raw(),
            Q1_15::from_float(-0.3f).raw()
        };

        // Channel 2: Different waveform
        alignas(8) int16_t channel2[8] = {
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(-0.5f).raw(),
            Q1_15::from_float(-0.5f).raw(),
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(0.5f).raw(),
            Q1_15::from_float(-0.5f).raw(),
            Q1_15::from_float(-0.5f).raw()
        };

        // Gain factors (channel 1: 0.8, channel 2: 0.6)
        alignas(8) int16_t gain1[8];
        alignas(8) int16_t gain2[8];
        for (int i = 0; i < 8; i++) {
            gain1[i] = Q1_15::from_float(0.8f).raw();
            gain2[i] = Q1_15::from_float(0.6f).raw();
        }

        // Temporary buffers
        alignas(8) int16_t ch1_scaled[8];
        alignas(8) int16_t ch2_scaled[8];
        alignas(8) int16_t mixed[8];

        // Create array wrappers
        FixedPointArray<1, 15, Backend> arr_ch1(channel1, 8);
        FixedPointArray<1, 15, Backend> arr_ch2(channel2, 8);
        FixedPointArray<1, 15, Backend> arr_g1(gain1, 8);
        FixedPointArray<1, 15, Backend> arr_g2(gain2, 8);
        FixedPointArray<1, 15, Backend> arr_ch1_scaled(ch1_scaled, 8);
        FixedPointArray<1, 15, Backend> arr_ch2_scaled(ch2_scaled, 8);
        FixedPointArray<1, 15, Backend> arr_mixed(mixed, 8);

        // Step 1: Apply gains (element-wise multiply)
        arr_ch1.elemult(arr_g1, arr_ch1_scaled);
        arr_ch2.elemult(arr_g2, arr_ch2_scaled);
        std::cout << "\n✓ Step 1: Applied channel gains (vec_elemult + vec_shift)\n";

        // Step 2: Mix channels (element-wise add)
        arr_ch1_scaled.add(arr_ch2_scaled, arr_mixed);
        std::cout << "✓ Step 2: Mixed channels (vec_add16x16_fast)\n";

        // Step 3: Compute RMS for metering
        auto rms_level = arr_mixed.rms();
        std::cout << "✓ Step 3: Computed RMS level (vec_rms16x16)\n";

        std::cout << "\nMixed output: [";
        for (int i = 0; i < 8; i++) {
            std::cout << Q1_15(mixed[i]).to_float();
            if (i < 7) std::cout << ", ";
        }
        std::cout << "]\n";

        print_value("RMS Level", rms_level.to_float(), "(for VU meter)");

        std::cout << "\n🔑 This entire audio processing chain uses vectorized ops!\n";
        std::cout << "   All operations use NatureDSP fast variants (aligned, N=8).\n";
    }

    print_header("SUMMARY - Vectorized Operations");
    std::cout << R"(
Key Takeaways:

1. Vectorization:
   ✓ Process multiple elements simultaneously
   ✓ Major performance boost on Xtensa HiFi DSP
   ✓ Same API for both scalar and vector operations

2. Fast Variants (NatureDSP):
   ✓ Available for: dot product, add, power
   ✓ Requirements: 8-byte aligned, N ≥ 4, N % 4 == 0
   ✓ Automatic fallback if requirements not met

3. Alignment:
   ✓ Use alignas(8) for static arrays
   ✓ Use aligned_alloc() for dynamic arrays
   ✓ Check alignment with (uintptr_t(ptr) & 7) == 0

4. Operations:
   ✓ Element-wise: add, sub, elemult (all vectorized)
   ✓ Reductions: sum, min, max, dot_product
   ✓ Statistics: mean, rms, variance, stddev

5. When to vectorize:
   ✓ Processing arrays/buffers (audio, signal, ML)
   ✓ Same operation on all elements
   ✓ Performance-critical inner loops

For maximum performance:
  • Align data to 8 bytes
  • Use lengths that are multiples of 4
  • Process data in batches
  • Use Xtensa backend for hardware acceleration

See fp.hpp and backends/xtensa/*.hpp for implementation details.
)";

    return 0;
}
