#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Trigonometric Implementation with Priority Dispatch
// ============================================================================
//
// Trigonometric functions use priority_tag dispatch:
//   - Priority 1: 16-bit operations (could use NatureDSP vec_sin16/cos16/etc.)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== SIN ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_sin_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template sin<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_sin_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_sin16 when available
    // For now, use portable implementation
    return xtensa_sin_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_sin_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_sin_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== COS ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_cos_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template cos<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_cos_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_cos16 when available
    // For now, use portable implementation
    return xtensa_cos_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_cos_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_cos_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== TAN ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_tan_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template tan<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_tan_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_tan16 when available
    // For now, use portable implementation
    return xtensa_tan_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_tan_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_tan_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== ATAN ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_atan_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template atan<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_atan_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_atan16 when available
    // For now, use portable implementation
    return xtensa_atan_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_atan_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_atan_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
