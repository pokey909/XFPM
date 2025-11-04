#pragma once
#include "../../helpers.hpp"
#include <cmath>

namespace fp {

// Forward declare the backend struct
struct ReferenceBackend;

// Power operation implementation for ReferenceBackend
namespace detail {

// Power function: base^exponent
// Base and exponent can have different Q formats
// Returns 0 for base <= 0
template<int Xb, int Yb>
inline typename StorageForBits<Xb>::type
reference_pow(typename StorageForBits<Xb>::type base,
              typename StorageForBits<Yb>::type exponent,
              int base_frac_bits,
              int exp_frac_bits)
{
    using Out = typename StorageForBits<Xb>::type;

    // Check for negative or zero base
    if (base <= 0) return 0;

    // Convert base to float
    float base_f = static_cast<float>(base) / static_cast<float>(1u << base_frac_bits);

    // Convert exponent to float
    float exp_f = static_cast<float>(exponent) / static_cast<float>(1u << exp_frac_bits);

    // Compute pow(base, exponent)
    float result = std::pow(base_f, exp_f);

    // Convert back to fixed-point with same format as base
    long long result_scaled = llroundf(result * static_cast<float>(1u << base_frac_bits));

    return sat_cast<Out>(result_scaled);
}

} // namespace detail
} // namespace fp
