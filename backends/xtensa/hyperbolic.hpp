#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Hyperbolic Implementation with Priority Dispatch
// ============================================================================
//
// Hyperbolic functions use priority_tag dispatch:
//   - Priority 1: 16-bit operations (could use NatureDSP vec_tanh16)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== TANH ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_tanh_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template tanh<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_tanh_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_tanh16 when available
    // For now, use portable implementation
    return xtensa_tanh_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_tanh_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_tanh_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
