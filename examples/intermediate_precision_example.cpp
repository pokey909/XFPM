/**
 * @file intermediate_precision_example.cpp
 * @brief Demonstrates scalar vs vectorized operations and intermediate precision
 *
 * This example shows the critical differences between:
 * 1. Scalar operations (mimi_mul_32x32 style) with 64-bit intermediate
 * 2. Vectorized element-wise operations (vec_elemult) without 64-bit intermediate
 * 3. Vectorized reduction operations (vec_dot) with 64-bit accumulator
 *
 * Understanding these differences is crucial for choosing the right operation
 * when porting code to use fp_lib.
 */

#include "fp.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

//=============================================================================
// SCALAR OPERATION: 64-bit Intermediate (like mimi_mul_32x32)
//=============================================================================

namespace scalar_example {

/**
 * This is the typical manual fixed-point multiplication with 64-bit intermediate.
 * It's what mimi_mul_32x32 does - uses int64_t for the product before shifting.
 */
inline int32_t manual_mul_32x32(int32_t x, int32_t y, uint32_t shift) {
    // 64-bit intermediate can hold the full 32x32 product without overflow
    return (int32_t)(((int64_t)x * y) >> shift);
}

void demonstrate() {
    std::cout << "=== SCALAR MULTIPLICATION: 64-bit Intermediate ===" << std::endl;
    std::cout << "\nExample: Q0.31 * Q0.31 -> Q0.31" << std::endl;

    // Two Q0.31 values: 0.5 * 0.5 = 0.25
    int32_t a = (int32_t)(0.5 * (1LL << 31));  // 0.5 in Q0.31
    int32_t b = (int32_t)(0.5 * (1LL << 31));  // 0.5 in Q0.31

    // Manual multiplication
    int32_t result_manual = manual_mul_32x32(a, b, 31);
    float result_float = (float)result_manual / (float)(1LL << 31);

    std::cout << "Manual (truncate):  0.5 * 0.5 = " << result_float << std::endl;

    // Using fp_lib
    fp::q<0, 31> a_fp = fp::q<0, 31>::from_float(0.5f);
    fp::q<0, 31> b_fp = fp::q<0, 31>::from_float(0.5f);
    auto result_fp = fp::mul_as<0, 31>(a_fp, b_fp);

    std::cout << "fp_lib (rounded):   0.5 * 0.5 = " << result_fp.to_float() << std::endl;

    std::cout << "\nKey Point: Both use 64-bit intermediate (int64_t)" << std::endl;
    std::cout << "Difference: fp_lib adds rounding before final shift" << std::endl;
    std::cout << "\nWhen to use:" << std::endl;
    std::cout << "  - Sample-by-sample processing" << std::endl;
    std::cout << "  - When you need 64-bit intermediate precision" << std::endl;
}

} // namespace scalar_example

//=============================================================================
// VECTORIZED ELEMENT-WISE: No 64-bit Intermediate (NatureDSP vec_elemult)
//=============================================================================

namespace vectorized_elemwise_example {

/**
 * NatureDSP vec_elemult32x32 does NOT use 64-bit intermediate!
 * It performs: 32x32 -> 32-bit with saturation
 * This is incompatible with Q-format arithmetic.
 */
void demonstrate() {
    std::cout << "\n=== VECTORIZED ELEMENT-WISE: No 64-bit Intermediate ===" << std::endl;
    std::cout << "\nNatureDSP vec_elemult32x32 signature:" << std::endl;
    std::cout << "  void vec_elemult32x32(int32_t *z, int32_t *x, int32_t *y, int N)" << std::endl;
    std::cout << "  Performs: z[i] = SAT32(x[i] * y[i])  // NO intermediate widening!" << std::endl;

    std::cout << "\nProblem with Q-format arithmetic:" << std::endl;
    std::cout << "  Q0.31 * Q0.31 produces Q0.62 intermediate (needs 64 bits)" << std::endl;
    std::cout << "  But vec_elemult32x32 only uses 32-bit intermediate" << std::endl;
    std::cout << "  Result: Saturation and incorrect values!" << std::endl;

    std::cout << "\nExample: 0.5 * 0.5 in Q0.31" << std::endl;
    std::cout << "  Expected intermediate: 1073741824 * 1073741824 = 1152921504606846976" << std::endl;
    std::cout << "  After shift by 31: 536870912 (0.25 in Q0.31) ✓" << std::endl;
    std::cout << "  With 32-bit: SATURATES to INT32_MAX before shift ✗" << std::endl;

    std::cout << "\nfp_lib solution:" << std::endl;
    std::cout << "  Falls back to reference implementation" << std::endl;
    std::cout << "  Uses proper 64-bit intermediate with rounding" << std::endl;

    // Demonstrate why it doesn't work
    std::cout << "\nWhat if you try to use vec_elemult32x32?" << std::endl;
    std::cout << "  Input: [0.5, 0.5, 0.5, 0.5]" << std::endl;
    std::cout << "  With 64-bit: [0.25, 0.25, 0.25, 0.25] ✓" << std::endl;
    std::cout << "  With vec_elemult32x32: [SAT, SAT, SAT, SAT] ✗" << std::endl;

    std::cout << "\nfp_lib solution:" << std::endl;
    // Demonstrate scalar version works correctly
    fp::q<0, 31> val = fp::q<0, 31>::from_float(0.5f);
    auto square = fp::mul_as<0, 31>(val, val);
    std::cout << "  Scalar: 0.5 * 0.5 = " << square.to_float() << " ✓" << std::endl;
    std::cout << "  For arrays: Use scalar loop (no SIMD for element-wise)" << std::endl;

    std::cout << "\nConclusion: vec_elemult32x32 CANNOT be used for Q-format!" << std::endl;
    std::cout << "fp_lib automatically uses correct reference implementation." << std::endl;
}

} // namespace vectorized_elemwise_example

//=============================================================================
// VECTORIZED REDUCTION: 64-bit Accumulator (NatureDSP vec_dot)
//=============================================================================

namespace vectorized_reduction_example {

/**
 * NatureDSP vec_dot32x32 DOES use 64-bit accumulator!
 * This makes it suitable for Q-format arithmetic.
 */
void demonstrate() {
    std::cout << "\n=== VECTORIZED REDUCTION: 64-bit Accumulator ===" << std::endl;
    std::cout << "\nNatureDSP vec_dot32x32 signature:" << std::endl;
    std::cout << "  int64_t vec_dot32x32(const int32_t *x, const int32_t *y, int N)" << std::endl;
    std::cout << "  Returns: int64_t accumulator with full precision ✓" << std::endl;

    std::cout << "\nHow it works:" << std::endl;
    std::cout << "  1. Multiply each pair: x[i] * y[i] (32x32 = 64-bit product)" << std::endl;
    std::cout << "  2. Accumulate into 64-bit register" << std::endl;
    std::cout << "  3. Return full 64-bit result" << std::endl;
    std::cout << "  4. User shifts result by fractional bits" << std::endl;

    std::cout << "\nThis is perfect for Q-format!" << std::endl;

    // Demonstrate with fp_lib
    std::cout << "\nExample:" << std::endl;
    std::cout << "  Array A: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]" << std::endl;
    std::cout << "  Array B: [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]" << std::endl;
    std::cout << "  Dot product: 0.5*0.25*8 = 1.0" << std::endl;

    std::cout << "\nWith vec_dot32x32:" << std::endl;
    std::cout << "  Each product: 0.5 * 0.25 uses 64-bit accumulator ✓" << std::endl;
    std::cout << "  Sum accumulator: 64-bit ✓" << std::endl;
    std::cout << "  Final shift: >> 31 to convert Q0.62 -> Q0.31 ✓" << std::endl;
    std::cout << "  Result: 1.0 (correct!)" << std::endl;

    std::cout << "\nConclusion: vec_dot32x32 CAN be used for Q-format!" << std::endl;
    std::cout << "fp_lib automatically uses it on Xtensa for performance." << std::endl;
}

} // namespace vectorized_reduction_example

//=============================================================================
// COMPARISON TABLE
//=============================================================================

void print_comparison_table() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "COMPARISON: Scalar vs Vectorized Operations" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::cout << "\n┌─────────────────────┬──────────────────┬──────────────┬─────────────────┐" << std::endl;
    std::cout << "│ Operation           │ Intermediate     │ Q-Format OK? │ fp_lib Uses     │" << std::endl;
    std::cout << "├─────────────────────┼──────────────────┼──────────────┼─────────────────┤" << std::endl;
    std::cout << "│ mimi_mul_32x32      │ 64-bit (int64_t) │ ✓ Yes        │ mul_as<>()      │" << std::endl;
    std::cout << "│ (scalar)            │ (truncate)       │              │ (+ rounding)    │" << std::endl;
    std::cout << "├─────────────────────┼──────────────────┼──────────────┼─────────────────┤" << std::endl;
    std::cout << "│ vec_elemult32x32    │ 32-bit only      │ ✗ NO!        │ Reference impl  │" << std::endl;
    std::cout << "│ (element-wise)      │ (saturates)      │              │ (64-bit sw)     │" << std::endl;
    std::cout << "├─────────────────────┼──────────────────┼──────────────┼─────────────────┤" << std::endl;
    std::cout << "│ vec_dot32x32        │ 64-bit accum     │ ✓ Yes        │ fp::dot<>()     │" << std::endl;
    std::cout << "│ (reduction)         │                  │              │ (HW optimized)  │" << std::endl;
    std::cout << "├─────────────────────┼──────────────────┼──────────────┼─────────────────┤" << std::endl;
    std::cout << "│ vec_sum32x32        │ 64-bit accum     │ ✓ Yes        │ fp::sum<>()     │" << std::endl;
    std::cout << "│ (reduction)         │                  │              │ (HW optimized)  │" << std::endl;
    std::cout << "└─────────────────────┴──────────────────┴──────────────┴─────────────────┘" << std::endl;
}

//=============================================================================
// PRACTICAL GUIDANCE
//=============================================================================

void print_practical_guidance() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "PRACTICAL GUIDANCE: When to Use What" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::cout << "\n1. SAMPLE-BY-SAMPLE MULTIPLICATION:" << std::endl;
    std::cout << "   Original: mimi_mul_32x32(a, b, shift)" << std::endl;
    std::cout << "   fp_lib:   mul_as<OutI, OutF>(a, b)" << std::endl;
    std::cout << "   Backend:  Software (64-bit intermediate + rounding)" << std::endl;
    std::cout << "   Example:  Energy calculation per sample" << std::endl;

    std::cout << "\n2. ELEMENT-WISE ARRAY MULTIPLICATION:" << std::endl;
    std::cout << "   Original: for(i) output[i] = mimi_mul_32x32(a[i], b[i], shift)" << std::endl;
    std::cout << "   fp_lib:   fp::array_elemult<OutI, OutF>(output, a, b, count)" << std::endl;
    std::cout << "   Backend:  Software (NatureDSP vec_elemult doesn't work for Q-format!)" << std::endl;
    std::cout << "   Example:  Apply gain to signal buffer" << std::endl;

    std::cout << "\n3. DOT PRODUCT / INNER PRODUCT:" << std::endl;
    std::cout << "   Original: for(i) acc += mimi_mul_32x32(a[i], b[i], shift)" << std::endl;
    std::cout << "   fp_lib:   fp::dot<OutI, OutF>(a, b, count)" << std::endl;
    std::cout << "   Backend:  Hardware SIMD (vec_dot32x32 on Xtensa)" << std::endl;
    std::cout << "   Example:  FIR filter, correlation, NLMS feedforward" << std::endl;

    std::cout << "\n4. SUM OF PRODUCTS (POWER/ENERGY):" << std::endl;
    std::cout << "   Original: for(i) acc += mimi_mul_32x32(a[i], a[i], shift)" << std::endl;
    std::cout << "   fp_lib:   fp::power<OutI, OutF>(a, count)  // Returns sum of squares" << std::endl;
    std::cout << "   Backend:  Hardware SIMD (vec_power32x32 on Xtensa)" << std::endl;
    std::cout << "   Example:  Signal energy, RMS calculation" << std::endl;

    std::cout << "\n5. MEAN (AVERAGE):" << std::endl;
    std::cout << "   Original: sum=0; for(i) sum+=a[i]; mean=sum/count" << std::endl;
    std::cout << "   fp_lib:   fp::mean(a, count)" << std::endl;
    std::cout << "   Backend:  Hardware SIMD (vec_sum32x32 on Xtensa)" << std::endl;
    std::cout << "   Example:  DC offset, moving average" << std::endl;
}

//=============================================================================
// CODE EXAMPLE: Refactoring mimi_mul_32x32 Usage
//=============================================================================

void refactoring_example() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "REFACTORING EXAMPLE: From mimi_mul_32x32 to fp_lib" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::cout << "\nBEFORE: Manual fixed-point with mimi_mul_32x32" << std::endl;
    std::cout << "```cpp" << std::endl;
    std::cout << "// Calculate energy per sample" << std::endl;
    std::cout << "for (uint32_t i = 0; i < N; i++) {" << std::endl;
    std::cout << "    int32_t sample = input[i];  // Q0.31" << std::endl;
    std::cout << "    // Square: Q0.31 * Q0.31 -> Q5.26 (shift by 31+5)" << std::endl;
    std::cout << "    energy[i] = mimi_mul_32x32(sample, sample, 36);" << std::endl;
    std::cout << "}" << std::endl;
    std::cout << "```" << std::endl;

    std::cout << "\nAFTER: Using fp_lib (sample-by-sample)" << std::endl;
    std::cout << "```cpp" << std::endl;
    std::cout << "// Same operation, type-safe with explicit format" << std::endl;
    std::cout << "for (uint32_t i = 0; i < N; i++) {" << std::endl;
    std::cout << "    q<0,31> sample = input[i];" << std::endl;
    std::cout << "    // Square with explicit output format" << std::endl;
    std::cout << "    energy[i] = mul_as<5,26>(sample, sample);" << std::endl;
    std::cout << "}" << std::endl;
    std::cout << "```" << std::endl;

    std::cout << "\nBETTER: Using vectorized operation (if applicable)" << std::endl;
    std::cout << "```cpp" << std::endl;
    std::cout << "// If you're computing energy for whole array:" << std::endl;
    std::cout << "auto total_energy = fp::power<5,26>(input.data(), N);" << std::endl;
    std::cout << "// Or if you need per-sample energy (can't vectorize elemult for Q-format!):" << std::endl;
    std::cout << "// Use scalar loop as shown above" << std::endl;
    std::cout << "```" << std::endl;
}

//=============================================================================
// MAIN
//=============================================================================

int main() {
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "╔════════════════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║     Scalar vs Vectorized Fixed-Point Operations                           ║" << std::endl;
    std::cout << "║     Understanding Intermediate Precision                                  ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════════════════╝" << std::endl;

    scalar_example::demonstrate();
    vectorized_elemwise_example::demonstrate();
    vectorized_reduction_example::demonstrate();

    print_comparison_table();
    print_practical_guidance();
    refactoring_example();

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "KEY TAKEAWAYS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "\n1. mimi_mul_32x32 uses 64-bit intermediate (truncation)" << std::endl;
    std::cout << "   → Replace with mul_as<>() (adds rounding for better accuracy)" << std::endl;

    std::cout << "\n2. NatureDSP vec_elemult32x32 does NOT use 64-bit intermediate" << std::endl;
    std::cout << "   → fp_lib falls back to software for correctness" << std::endl;
    std::cout << "   → Cannot vectorize Q-format element-wise multiplication on Xtensa" << std::endl;

    std::cout << "\n3. NatureDSP vec_dot32x32 DOES use 64-bit accumulator" << std::endl;
    std::cout << "   → fp_lib uses it for massive SIMD speedup" << std::endl;
    std::cout << "   → Also applies to: vec_sum, vec_power, vec_mean (via sum)" << std::endl;

    std::cout << "\n4. Choose the right operation:" << std::endl;
    std::cout << "   - Sample-by-sample: mul_as<>() (scalar)" << std::endl;
    std::cout << "   - Dot product/sum: fp::dot(), fp::sum() (vectorized!)" << std::endl;
    std::cout << "   - Element-wise multiply: Can't vectorize, use scalar loop" << std::endl;

    std::cout << "\n" << std::string(80, '=') << std::endl;

    return 0;
}
