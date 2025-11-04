#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Power Implementation with Priority Dispatch
// ============================================================================
//
// This implementation uses priority_tag dispatch to select specialized
// power paths based on operand bucket sizes:
//   - Priority 1: 32-bit base (uses NatureDSP vec_pow_32x32 when available)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order (Priority 0 first,
// then Priority 1) so that each level can call the next without forward declarations.

// -------- Priority 0: Generic Fallback → ReferenceBackend (define first) --------

template<int Xb, int Yb>
inline Storage_t<Xb>
xtensa_pow_impl(Storage_t<Xb> base,
                Storage_t<Yb> exponent,
                int base_frac_bits,
                int exp_frac_bits,
                priority_tag<0>)
{
    // Fallback to portable ReferenceBackend implementation
    return ReferenceBackend::template pow<Xb, Yb>(base, exponent, base_frac_bits, exp_frac_bits);
}

// -------- Priority 1: 32-bit base Specialization --------

template<int Xb>
using is_base_32 = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

// Enabled when base is 32-bit
template<int Xb, int Yb, EnableIf<is_base_32<Xb>::value> = 0>
inline int32_t
xtensa_pow_impl(int32_t base,
                Storage_t<Yb> exponent,
                int base_frac_bits,
                int exp_frac_bits,
                priority_tag<1>)
{
    // TODO: Use NatureDSP vec_pow_32x32 function when available:
    // For now, use portable implementation
    return ReferenceBackend::template pow<Xb, Yb>(base, exponent, base_frac_bits, exp_frac_bits);
}

// Forward to Priority 0 when base is NOT 32-bit
template<int Xb, int Yb, EnableIf<!is_base_32<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_pow_impl(Storage_t<Xb> base,
                Storage_t<Yb> exponent,
                int base_frac_bits,
                int exp_frac_bits,
                priority_tag<1>)
{
    return xtensa_pow_impl<Xb, Yb>(base, exponent, base_frac_bits, exp_frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
