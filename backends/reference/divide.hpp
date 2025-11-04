#pragma once
#include "../../helpers.hpp"

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Division Implementation
// ============================================================================
//
// Division in fixed-point:
//   If a is Qx.A and b is Qy.B, and we want output Qz.C:
//   - Raw division: (a / b) has implicit fractional bits = A - B
//   - To get C fractional bits, we need to shift the dividend:
//     result = (a << shift) / b, where shift = C - (A - B) = C - A + B
//
// For our implementation:
//   - Dividend 'ax' has F_a fractional bits
//   - Divisor 'by' has F_b fractional bits
//   - Output should have F_out fractional bits
//   - shift = F_out - (F_a - F_b) = F_out - F_a + F_b
//
// However, fp.hpp passes us 'out_frac_shift' which is calculated differently.
// Looking at multiply: out_frac_shift = frac_in - frac_out = (F_a + F_b) - F_out
// For divide, we need the opposite operation, so:
//   If we follow the same pattern, out_frac_shift = (F_a - F_b) - F_out
//   Then actual_shift = -out_frac_shift gives us what we need.

template<int Xb, int Yb, int Ob>
inline typename StorageForBits<Ob>::type
reference_div( typename StorageForBits<Xb>::type ax,
               typename StorageForBits<Yb>::type by,
               int out_frac_shift )
{
    using Out = typename StorageForBits<Ob>::type;

    // Handle division by zero
    if (by == 0) {
        // Return saturated max/min based on dividend sign
        if (ax >= 0) {
            return std::numeric_limits<Out>::max();
        } else {
            return std::numeric_limits<Out>::min();
        }
    }

    // For division, we need to shift left before dividing
    // The shift amount is the negative of what multiply uses
    int shift = -out_frac_shift;

    // Use wide intermediate to avoid overflow
    long long dividend = static_cast<long long>(ax);
    long long divisor = static_cast<long long>(by);

    // Apply shift to dividend
    if (shift >= 0) {
        dividend = dividend << shift;
    } else {
        // If shift is negative, we need to shift the divisor instead
        divisor = divisor << (-shift);
    }

    // Perform division with rounding
    // For rounding towards nearest, add half of divisor before division
    long long quotient;
    if ((dividend >= 0 && divisor > 0) || (dividend < 0 && divisor < 0)) {
        // Same sign: round towards +infinity
        quotient = (dividend + divisor/2) / divisor;
    } else {
        // Opposite sign: round towards -infinity
        quotient = (dividend - divisor/2) / divisor;
    }

    return sat_cast<Out>(quotient);
}

} // namespace detail
} // namespace fp
