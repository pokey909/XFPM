#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Logarithm Implementation with Priority Dispatch
// ============================================================================
//
// Logarithm functions use NatureDSP library functions which work with 32-bit
// Q16.15 input and Q6.25 output:
//   - Priority 1: 32-bit operations (uses NatureDSP scl_log2_32x32, etc.)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== LOG2 ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_log2_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template log2<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb,32>::value>;

// Enabled when 32-bit
template<int Xb, typename std::enable_if<is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_log2_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_log2_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_log2_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, typename std::enable_if<!is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_log2_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<1>)
{
    return xtensa_log2_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== LOGN (Natural Log) ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_logn_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template logn<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

// Enabled when 32-bit
template<int Xb, typename std::enable_if<is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_logn_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_logn_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_logn_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, typename std::enable_if<!is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_logn_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<1>)
{
    return xtensa_logn_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== LOG10 ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline int32_t
xtensa_log10_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template log10<Xb>(ax, frac_bits);
}

// -------- Priority 1: 32-bit Specialization --------

// Enabled when 32-bit
template<int Xb, typename std::enable_if<is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_log10_impl(int32_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_log10_32x32 when available
    // For now, fallback to portable implementation
    return xtensa_log10_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 32-bit
template<int Xb, typename std::enable_if<!is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_log10_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<1>)
{
    return xtensa_log10_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
