#pragma once
#include "../../helpers.hpp"
#include <cmath>
#include <limits>

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Square Root Implementation
// ============================================================================
//
// Square root functions follow NatureDSP conventions:
//   - Input and output have the same Q format
//   - Returns 0x80000000/0x8000 (most negative value) on negative input
//
// Implementation:
//   - For Q<I,F> format: sqrt(x_raw / 2^F) = result / 2^F
//   - This means: sqrt(x_raw) = result * sqrt(2^F) = result * 2^(F/2)
//   - We compute: result = sqrt(x_raw) / 2^(F/2)
//
// Algorithm:
//   1. Check for negative input (return error value)
//   2. Convert fixed-point to float
//   3. Compute floating-point sqrt
//   4. Convert back to same Q format with proper scaling

template<int Xb>
inline typename StorageForBits<Xb>::type
reference_sqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
{
    using storage_t = typename StorageForBits<Xb>::type;

    // Check for negative input
    if (ax < 0) {
        return std::numeric_limits<storage_t>::min();  // 0x80000000 or 0x8000
    }

    // Convert to float: value = ax / 2^frac_bits
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute sqrt(x)
    float result = std::sqrt(x);

    // Convert result back to Q format: result_raw = result * 2^frac_bits
    float scale = static_cast<float>(1u << frac_bits);
    long long result_scaled = llroundf(result * scale);

    return sat_cast<storage_t>(result_scaled);
}

} // namespace detail
} // namespace fp
