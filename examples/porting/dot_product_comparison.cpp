/**
 * @file dot_product_comparison.cpp
 * @brief Side-by-side comparison of manual fixed-point vs fp_lib
 *
 * This example shows a real NLMS dot product operation implemented both ways.
 * The manual version is adapted from nlms_node.cpp:320-327
 */

#include <cstdint>
#include <array>
#include <iostream>
#include <cmath>

// For the refactored version
#include "fp.hpp"

//=============================================================================
// BEFORE: Manual Fixed-Point Implementation
//=============================================================================

namespace manual {

/**
 * Helper function for Q31 fixed-point multiplication
 * Multiplies two Q-format numbers and returns result in specified Q-format
 * NOTE: This is the ACTUAL implementation - NO rounding, just truncation!
 */
inline int32_t mimi_mul_32x32(int32_t a, int32_t b, int fractional_bits) {
    // trust the caller that the result won't overflow
    return (int32_t)(((int64_t)a * b) >> fractional_bits);
}

/**
 * Dot product: buffer_a (Q31) * buffer_b (Q26) -> Q26
 *
 * ISSUES:
 * - Q-formats only documented in comments
 * - Caller must know buffer_a is Q31, buffer_b is Q26
 * - Magic number 31 for fractional bits
 * - No compile-time type safety
 * - Overflow/saturation handling is manual
 */
int32_t dot_product(const int32_t *buffer_a, const int32_t *buffer_b, uint32_t count) {
    int32_t accumulator = 0;
    for (uint32_t idx = 0; idx < count; idx++) {
        // Q31 * Q26 with 31 fractional bits -> Q26
        accumulator += mimi_mul_32x32(buffer_a[idx], buffer_b[idx], 31);
    }
    return accumulator;  // Result is Q26
}

// Example usage
void example_manual() {
    constexpr uint32_t LENGTH = 8;

    // Speaker buffer in Q31 format (comment documentation only!)
    int32_t speaker_buffer[LENGTH];

    // IR estimate in Q26 format (comment documentation only!)
    int32_t ir_estimate[LENGTH];

    // Initialize with some values
    for (uint32_t i = 0; i < LENGTH; i++) {
        // 0.5 in Q31: 0.5 * 2^31 = 1073741824
        speaker_buffer[i] = 1073741824;
        // 0.25 in Q26: 0.25 * 2^26 = 16777216
        ir_estimate[i] = 16777216;
    }

    // Call dot product - hope we got the formats right!
    int32_t result = dot_product(speaker_buffer, ir_estimate, LENGTH);

    // Convert Q26 back to float for display
    float result_float = (float)result / (float)(1 << 26);

    std::cout << "Manual implementation result: " << result_float << std::endl;
    std::cout << "  Expected: 0.5 * 0.25 * 8 = 1.0" << std::endl;
}

} // namespace manual

//=============================================================================
// AFTER: Using fp_lib
//=============================================================================

namespace with_fplib {

/**
 * The library already provides dot product as a built-in function!
 *
 * BENEFITS:
 * - Built-in function - no need to implement yourself!
 * - Q-formats are explicit in type signatures
 * - Compile-time verification of format compatibility
 * - Automatic saturation and rounding
 * - No magic numbers
 * - Self-documenting code
 * - Backend-optimized (uses SIMD when available)
 */

// Example usage with explicit types matching NLMS
void example_fplib() {
    constexpr uint32_t LENGTH = 8;

    // Speaker buffer in Q0.31 format (explicit in type, matches original Q31)
    std::array<q<0, 31>, LENGTH> speaker_buffer;

    // IR estimate in Q5.26 format (explicit in type, matches original Q26)
    std::array<q<5, 26>, LENGTH> ir_estimate;

    // Initialize with type-safe from_float()
    for (uint32_t i = 0; i < LENGTH; i++) {
        speaker_buffer[i] = q<0, 31>::from_float(0.5f);
        ir_estimate[i] = q<5, 26>::from_float(0.25f);
    }

    // Use built-in dot product with explicit output format
    // Q0.31 * Q5.26 -> Q5.26 (matches original behavior, stays in 32 bits)
    auto result = fp::dot<5, 26>(speaker_buffer.data(), ir_estimate.data(), LENGTH);

    // Type-safe conversion back to float
    float result_float = result.to_float();

    std::cout << "fp_lib implementation result: " << result_float << std::endl;
    std::cout << "  Result type: q<" << result.I << "," << result.F << "> (32-bit format)" << std::endl;
    std::cout << "  Expected: 0.5 * 0.25 * 8 = 1.0" << std::endl;
}

} // namespace with_fplib

//=============================================================================
// COMPARISON: Energy Calculation
//=============================================================================

namespace energy_comparison {

// BEFORE: Manual fixed-point energy calculation
namespace manual {

inline int32_t mimi_mul_32x32(int32_t a, int32_t b, int fractional_bits) {
    // trust the caller that the result won't overflow
    return (int32_t)(((int64_t)a * b) >> fractional_bits);
}

void calculate_energy(int32_t mic_input, int32_t error,
                     int32_t& input_energy, int32_t& error_energy)
{
    // mic_input is Q31, error is Q26
    // Want output in Q26 format

    // Manual calculation with magic shift values
    input_energy = mimi_mul_32x32(mic_input, mic_input, 31 + 5);  // Q31 * Q31 -> Q26
    error_energy = mimi_mul_32x32(error, error, 26);               // Q26 * Q26 -> Q26

    // What if we need to change output format? Need to recalculate all shifts!
}

} // namespace manual

// AFTER: Using fp_lib with explicit format control
namespace with_fplib {

// Note: Original code stores energy in Q26 format (Q5.26), not higher precision
void calculate_energy(q<0, 31> mic_input, q<5, 26> error,
                     q<5, 26>& input_energy, q<5, 26>& error_energy)
{
    // Explicit output format matching original Q26 (Q5.26)
    // Q0.31 * Q0.31 -> Q5.26 (original shift: 31+5)
    input_energy = mul_as<5, 26>(mic_input, mic_input);
    // Q5.26 * Q5.26 -> Q5.26 (original shift: 26)
    error_energy = mul_as<5, 26>(error, error);

    // Output stays in 32-bit format, works on Xtensa!
}

} // namespace with_fplib

void demonstrate_energy() {
    std::cout << "\n=== Energy Calculation Comparison ===" << std::endl;

    // Test value: 0.5
    int32_t mic_q31 = (int32_t)(0.5f * (1 << 31));
    int32_t error_q26 = (int32_t)(0.5f * (1 << 26));

    // Manual version
    int32_t input_energy_manual, error_energy_manual;
    manual::calculate_energy(mic_q31, error_q26, input_energy_manual, error_energy_manual);

    float input_e_float = (float)input_energy_manual / (float)(1 << 26);
    std::cout << "Manual: input_energy = " << input_e_float << " (expected 0.25)" << std::endl;

    // fp_lib version
    q<5, 26> mic_fp = q<5, 26>::from_float(0.5f);
    q<5, 26> error_fp = q<5, 26>::from_float(0.5f);
    q<10, 52> input_energy_fp, error_energy_fp;

    with_fplib::calculate_energy(mic_fp, error_fp, input_energy_fp, error_energy_fp);

    std::cout << "fp_lib: input_energy = " << input_energy_fp.to_float() << " (expected 0.25)" << std::endl;
    std::cout << "  Output type: q<10,52> (automatically handled)" << std::endl;
}

} // namespace energy_comparison

//=============================================================================
// Main: Run all comparisons
//=============================================================================

int main() {
    std::cout << "=== Fixed-Point Library Refactoring Demo ===" << std::endl;
    std::cout << "\nThis demonstrates refactoring manual fixed-point code" << std::endl;
    std::cout << "from the NLMS algorithm to use the fp_lib library.\n" << std::endl;

    std::cout << "=== Dot Product Comparison ===" << std::endl;
    manual::example_manual();
    std::cout << std::endl;
    with_fplib::example_fplib();

    energy_comparison::demonstrate_energy();

    std::cout << "\n=== Key Advantages of fp_lib ===" << std::endl;
    std::cout << "✓ Type safety: Q-formats are part of the type system" << std::endl;
    std::cout << "✓ No magic numbers: No manual shift calculations" << std::endl;
    std::cout << "✓ Automatic saturation and rounding" << std::endl;
    std::cout << "✓ Self-documenting: q<5,26> vs int32_t with comment" << std::endl;
    std::cout << "✓ Easier refactoring: Change type declaration, not all shifts" << std::endl;
    std::cout << "✓ Backend system: Automatic platform optimization" << std::endl;

    return 0;
}
