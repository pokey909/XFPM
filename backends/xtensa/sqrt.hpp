#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Square Root Implementation with Priority Dispatch
// ============================================================================
//
// Square root functions use NatureDSP library functions:
//   - Priority 2: 32-bit operations (uses NatureDSP scl_sqrt32x32, etc.)
//   - Priority 1: 16-bit operations (uses NatureDSP scl_sqrt16x16, etc.)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline typename StorageForBits<Xb>::type
xtensa_sqrt_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template sqrt<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb,16>::value>;

// Enabled when 16-bit
template<int Xb, typename std::enable_if<is_16bit<Xb>::value, int>::type = 0>
inline int16_t
xtensa_sqrt_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_sqrt16x16 when available
    // For now, fallback to portable implementation
    return xtensa_sqrt_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, typename std::enable_if<!is_16bit<Xb>::value, int>::type = 0>
inline typename StorageForBits<Xb>::type
xtensa_sqrt_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<1>)
{
    return xtensa_sqrt_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// -------- Priority 2: 32-bit Specialization --------

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb,32>::value>;

// Enabled when 32-bit
template<int Xb, typename std::enable_if<is_32bit<Xb>::value, int>::type = 0>
inline int32_t
xtensa_sqrt_impl(int32_t ax, int frac_bits, priority_tag<2>)
{
    // TODO: Use NatureDSP scl_sqrt32x32 when available
    // For now, fallback to portable implementation
    return xtensa_sqrt_impl<Xb>(ax, frac_bits, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 32-bit
template<int Xb, typename std::enable_if<!is_32bit<Xb>::value, int>::type = 0>
inline typename StorageForBits<Xb>::type
xtensa_sqrt_impl(typename StorageForBits<Xb>::type ax, int frac_bits, priority_tag<2>)
{
    return xtensa_sqrt_impl<Xb>(ax, frac_bits, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
