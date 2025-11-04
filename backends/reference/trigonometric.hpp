#pragma once
#include "../../helpers.hpp"
#include <cmath>

namespace fp {
namespace detail {

// ============================================================================
// Reference Trigonometric Implementation
// ============================================================================
//
// Trigonometric functions using standard C++ math library.
// Input is treated as radians in the input Q format.
// Output maintains the same Q format as input.

// Sine function
template<int Xb>
inline Storage_t<Xb>
reference_sin(Storage_t<Xb> ax, int frac_bits)
{
    // Convert to float (radians)
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute sine
    float result = std::sin(x);

    // Convert back to fixed-point
    long long result_scaled = llroundf(result * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

// Cosine function
template<int Xb>
inline Storage_t<Xb>
reference_cos(Storage_t<Xb> ax, int frac_bits)
{
    // Convert to float (radians)
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute cosine
    float result = std::cos(x);

    // Convert back to fixed-point
    long long result_scaled = llroundf(result * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

// Tangent function
template<int Xb>
inline Storage_t<Xb>
reference_tan(Storage_t<Xb> ax, int frac_bits)
{
    // Convert to float (radians)
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute tangent
    float result = std::tan(x);

    // Convert back to fixed-point
    long long result_scaled = llroundf(result * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

// Arctangent function
template<int Xb>
inline Storage_t<Xb>
reference_atan(Storage_t<Xb> ax, int frac_bits)
{
    // Convert to float
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute arctangent (result in radians)
    float result = std::atan(x);

    // Convert back to fixed-point
    long long result_scaled = llroundf(result * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

} // namespace detail
} // namespace fp
