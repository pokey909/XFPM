#pragma once
#include "../../helpers.hpp"
#include <cmath>
#include <limits>

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Reciprocal Square Root Implementation
// ============================================================================
//
// Reciprocal square root computes 1/sqrt(x) for a fixed-point value.
//
// Unlike logarithm operations which convert to/from specific Q formats,
// rsqrt maintains the same Q format for input and output:
//   - If input is Qx.F, output is also Qx.F
//   - Returns error value (numeric_limits::min()) on negative or zero input
//
// Implementation strategy:
//   1. Convert input to floating-point
//   2. Compute 1/sqrt(x) using standard library
//   3. Convert result back to fixed-point with saturation
//
// Note: This follows the NatureDSP convention where rsqrt functions
// return the result in the same format as the input, with 1 LSB mantissa accuracy.

template<int Xb>
inline typename StorageForBits<Xb>::type
reference_rsqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
{
    using Out = typename StorageForBits<Xb>::type;

    // Check for negative or zero input - undefined behavior
    if (ax <= 0) {
        return std::numeric_limits<Out>::min();  // Return error value
    }

    // Convert from fixed-point to float
    // The value is ax / 2^frac_bits
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute reciprocal square root: 1/sqrt(x)
    float result = 1.0f / std::sqrt(x);

    // Convert back to fixed-point with rounding
    // Scale by 2^frac_bits and round to nearest
    float scale = static_cast<float>(1u << frac_bits);
    long long result_scaled = llroundf(result * scale);

    // Saturate to output type range
    return sat_cast<Out>(result_scaled);
}

} // namespace detail
} // namespace fp
