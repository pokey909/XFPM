#pragma once
#include "../../helpers.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Reference Vector Scale Implementation
// ============================================================================
//
// Array operations that modify elements in-place by shifting or scaling.
// These operations work on arrays of fixed-point values stored as integers.

// Shift all array elements by a bit amount (in-place)
// Positive shift_amount: left shift (multiply by 2^shift_amount)
// Negative shift_amount: right shift (divide by 2^|shift_amount|)
template<int Xb>
inline void
reference_array_shift(Storage_t<Xb>* arr, size_t length, int shift_amount)
{
    if (shift_amount == 0) {
        return;  // No shift needed
    }

    for (size_t i = 0; i < length; ++i) {
        long long val = static_cast<long long>(arr[i]);
        long long shifted;

        if (shift_amount > 0) {
            // Left shift
            shifted = val << shift_amount;
        } else {
            // Right shift (arithmetic)
            shifted = val >> (-shift_amount);
        }

        arr[i] = sat_cast<Storage_t<Xb>>(shifted);
    }
}

// Scale all array elements by a scalar fixed-point value (in-place)
// arr[i] = sat_cast((arr[i] * scale_factor) >> scale_frac_bits)
template<int Xb>
inline void
reference_array_scale(Storage_t<Xb>* arr, size_t length,
                      Storage_t<Xb> scale_factor, int scale_frac_bits)
{
    for (size_t i = 0; i < length; ++i) {
        long long val = static_cast<long long>(arr[i]);
        long long scale = static_cast<long long>(scale_factor);
        long long product = val * scale;
        long long scaled = round_shift(product, scale_frac_bits);
        arr[i] = sat_cast<Storage_t<Xb>>(scaled);
    }
}

} // namespace detail
} // namespace fp
