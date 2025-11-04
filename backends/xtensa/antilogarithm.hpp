#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Antilogarithm Implementation with Priority Dispatch
// ============================================================================
//
// Antilogarithm functions use NatureDSP library functions which work with
// Q6.25 input and Q16.15 output (OPPOSITE of logarithm):
//   - Priority 1: 32-bit operations (uses NatureDSP scl_antilog2_32x32, etc.)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== ANTILOG2 (2^x) ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_antilog2_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template antilog2<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb,32>::value>;

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilog2_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_antilog2_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_antilog2_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilog2_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_antilog2_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== ANTILOGN (e^x) ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_antilogn_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template antilogn<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilogn_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_antilogn_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_antilogn_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilogn_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_antilogn_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== ANTILOG10 (10^x) ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_antilog10_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template antilog10<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilog10_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_antilog10_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_antilog10_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_antilog10_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_antilog10_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
