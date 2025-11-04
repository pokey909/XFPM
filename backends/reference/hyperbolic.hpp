#pragma once
#include "../../helpers.hpp"
#include <cmath>

namespace fp {
namespace detail {

// ============================================================================
// Reference Hyperbolic Implementation
// ============================================================================
//
// Hyperbolic functions using standard C++ math library.
// Input is treated as a value in the input Q format.
// Output maintains the same Q format as input.

// Hyperbolic tangent function
template<int Xb>
inline Storage_t<Xb>
reference_tanh(Storage_t<Xb> ax, int frac_bits)
{
    // Convert to float
    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);

    // Compute hyperbolic tangent
    float result = std::tanh(x);

    // Convert back to fixed-point
    long long result_scaled = llroundf(result * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

} // namespace detail
} // namespace fp
