#pragma once
#include "../../helpers.hpp"
#include <cmath>
#include <limits>

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Logarithm Implementation
// ============================================================================
//
// Logarithm functions follow NatureDSP conventions:
//   - Input is interpreted as Q16.15 (regardless of actual format)
//   - Output is in Q6.25 format
//   - Returns 0x80000000 (most negative int32_t) on negative or zero input
//
// This means:
//   - Input value range: approximately [-1.0, 1.0)
//   - Output value range: approximately [-32, 32) with high precision
//
// Implementation uses floating-point logarithm functions and converts
// the result to Q6.25 fixed-point format.

// Base-2 logarithm (log2)
// Input: Any Q format, output: Q6.25
// The input is first converted to Q16.15 format as per NatureDSP convention
template<int Xb>
inline int32_t
reference_log2(typename StorageForBits<Xb>::type ax, int frac_bits)
{
    // First convert input to Q16.15 format
    // If input has more fractional bits than 15, shift right
    // If input has fewer fractional bits than 15, shift left
    int shift = frac_bits - 15;

    int32_t input_q16_15;
    if (shift > 0) {
        // More fractional bits than Q16.15, shift right
        input_q16_15 = static_cast<int32_t>(ax) >> shift;
    } else if (shift < 0) {
        // Fewer fractional bits than Q16.15, shift left
        input_q16_15 = static_cast<int32_t>(ax) << (-shift);
    } else {
        input_q16_15 = static_cast<int32_t>(ax);
    }

    // Check for negative or zero input
    if (input_q16_15 <= 0) {
        return std::numeric_limits<int32_t>::min();  // 0x80000000
    }

    // Convert Q16.15 to float
    float x = static_cast<float>(input_q16_15) / static_cast<float>(1u << 15);

    // Compute log2(x)
    float result = std::log2(x);

    // Convert result to Q6.25
    float scale_q6_25 = static_cast<float>(1u << 25);
    long long result_scaled = llroundf(result * scale_q6_25);

    return sat_cast<int32_t>(result_scaled);
}

// Natural logarithm (logn / ln)
template<int Xb>
inline int32_t
reference_logn(typename StorageForBits<Xb>::type ax, int frac_bits)
{
    int shift = frac_bits - 15;

    int32_t input_q16_15;
    if (shift > 0) {
        input_q16_15 = static_cast<int32_t>(ax) >> shift;
    } else if (shift < 0) {
        input_q16_15 = static_cast<int32_t>(ax) << (-shift);
    } else {
        input_q16_15 = static_cast<int32_t>(ax);
    }

    if (input_q16_15 <= 0) {
        return std::numeric_limits<int32_t>::min();
    }

    float x = static_cast<float>(input_q16_15) / static_cast<float>(1u << 15);
    float result = std::log(x);  // Natural logarithm

    float scale_q6_25 = static_cast<float>(1u << 25);
    long long result_scaled = llroundf(result * scale_q6_25);

    return sat_cast<int32_t>(result_scaled);
}

// Base-10 logarithm (log10)
template<int Xb>
inline int32_t
reference_log10(typename StorageForBits<Xb>::type ax, int frac_bits)
{
    int shift = frac_bits - 15;

    int32_t input_q16_15;
    if (shift > 0) {
        input_q16_15 = static_cast<int32_t>(ax) >> shift;
    } else if (shift < 0) {
        input_q16_15 = static_cast<int32_t>(ax) << (-shift);
    } else {
        input_q16_15 = static_cast<int32_t>(ax);
    }

    if (input_q16_15 <= 0) {
        return std::numeric_limits<int32_t>::min();
    }

    float x = static_cast<float>(input_q16_15) / static_cast<float>(1u << 15);
    float result = std::log10(x);

    float scale_q6_25 = static_cast<float>(1u << 25);
    long long result_scaled = llroundf(result * scale_q6_25);

    return sat_cast<int32_t>(result_scaled);
}

} // namespace detail
} // namespace fp
