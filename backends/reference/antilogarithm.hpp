#pragma once
#include "../../helpers.hpp"
#include <cmath>
#include <limits>

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Antilogarithm Implementation
// ============================================================================
//
// Antilogarithm functions follow NatureDSP conventions (OPPOSITE of logarithm):
//   - Input is interpreted as Q6.25 (regardless of actual format)
//   - Output is in Q16.15 format
//   - Returns 0x7FFFFFFF (max int32_t) on overflow
//   - Returns 0 on underflow
//
// This means:
//   - Input value range: approximately [-32, 32) with high precision
//   - Output value range: approximately [-32768, 32768)
//
// Implementation uses floating-point exponential functions and converts
// the result to Q16.15 fixed-point format.

// Base-2 antilogarithm (2^x)
// Input: Any Q format (interpreted as Q6.25), output: Q16.15
template<int Xb>
inline int32_t
reference_antilog2(Storage_t<Xb> ax, int frac_bits)
{
    // Input is interpreted as Q6.25 regardless of actual format
    // Convert input from Q6.25 to float
    float x = static_cast<float>(ax) / static_cast<float>(1u << 25);

    // Compute 2^x
    float result = std::pow(2.0f, x);

    // Convert result to Q16.15
    float scale_q16_15 = static_cast<float>(1u << 15);
    long long result_scaled = llroundf(result * scale_q16_15);

    // Handle overflow/underflow
    if (result_scaled > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<int32_t>::max();  // 0x7FFFFFFF
    }
    if (result_scaled < 0) {
        return 0;  // Underflow
    }

    return sat_cast<int32_t>(result_scaled);
}

// Natural antilogarithm (e^x)
// Input: Any Q format (interpreted as Q6.25), output: Q16.15
template<int Xb>
inline int32_t
reference_antilogn(Storage_t<Xb> ax, int frac_bits)
{
    // Input is interpreted as Q6.25 regardless of actual format
    // Convert input from Q6.25 to float
    float x = static_cast<float>(ax) / static_cast<float>(1u << 25);

    // Compute e^x
    float result = std::exp(x);

    // Convert result to Q16.15
    float scale_q16_15 = static_cast<float>(1u << 15);
    long long result_scaled = llroundf(result * scale_q16_15);

    // Handle overflow/underflow
    if (result_scaled > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<int32_t>::max();  // 0x7FFFFFFF
    }
    if (result_scaled < 0) {
        return 0;  // Underflow
    }

    return sat_cast<int32_t>(result_scaled);
}

// Base-10 antilogarithm (10^x)
// Input: Any Q format (interpreted as Q6.25), output: Q16.15
template<int Xb>
inline int32_t
reference_antilog10(Storage_t<Xb> ax, int frac_bits)
{
    // Input is interpreted as Q6.25 regardless of actual format
    // Convert input from Q6.25 to float
    float x = static_cast<float>(ax) / static_cast<float>(1u << 25);

    // Compute 10^x
    float result = std::pow(10.0f, x);

    // Convert result to Q16.15
    float scale_q16_15 = static_cast<float>(1u << 15);
    long long result_scaled = llroundf(result * scale_q16_15);

    // Handle overflow/underflow
    if (result_scaled > std::numeric_limits<int32_t>::max()) {
        return std::numeric_limits<int32_t>::max();  // 0x7FFFFFFF
    }
    if (result_scaled < 0) {
        return 0;  // Underflow
    }

    return sat_cast<int32_t>(result_scaled);
}

} // namespace detail
} // namespace fp
