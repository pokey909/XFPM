#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Reciprocal Square Root Implementation with Priority Dispatch
// ============================================================================
//
// This implementation uses priority_tag dispatch to select specialized
// rsqrt paths based on operand bucket size:
//   - Priority 2: 32-bit (uses NatureDSP scl_rsqrt32x32)
//   - Priority 1: 16-bit (uses NatureDSP scl_rsqrt16x16)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order (Priority 0 first,
// then Priority 1, then Priority 2) so that each level can call the next without forward declarations.
//
// NatureDSP variants:
//   - scl_rsqrt32x32: 32-bit reciprocal square root with 1 LSB mantissa accuracy
//   - scl_rsqrt16x16: 16-bit reciprocal square root with 1 LSB mantissa accuracy
//
// Input and output are in the same Q format.

// -------- Priority 0: Generic Fallback → ReferenceBackend (define first) --------

template<int Xb>
inline typename StorageForBits<Xb>::type
xtensa_rsqrt_impl(typename StorageForBits<Xb>::type ax,
                  int frac_bits,
                  priority_tag<0>)
{
    // Fallback to portable ReferenceBackend implementation
    return ReferenceBackend::template rsqrt<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_rsqrt_16 = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb,
         typename std::enable_if<is_rsqrt_16<Xb>::value, int>::type = 0>
inline int16_t
xtensa_rsqrt_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_rsqrt16x16 when available
    // For now, use portable implementation

    if (ax <= 0) {
        return std::numeric_limits<int16_t>::min();
    }

    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);
    float result = 1.0f / std::sqrt(x);
    float scale = static_cast<float>(1u << frac_bits);
    long long result_scaled = llroundf(result * scale);

    return sat_cast<int16_t>(result_scaled);
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb,
         typename std::enable_if<!is_rsqrt_16<Xb>::value, int>::type = 0>
inline typename StorageForBits<Xb>::type
xtensa_rsqrt_impl(typename StorageForBits<Xb>::type ax,
                  int frac_bits,
                  priority_tag<1>)
{
    return xtensa_rsqrt_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// -------- Priority 2: 32-bit Specialization --------

template<int Xb>
using is_rsqrt_32 = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

// Enabled when 32-bit
template<int Xb,
         typename std::enable_if<is_rsqrt_32<Xb>::value, int>::type = 0>
inline int32_t
xtensa_rsqrt_impl(int32_t ax, int frac_bits, priority_tag<2>)
{
    // TODO: Use NatureDSP scl_rsqrt32x32 when available
    // For now, use portable implementation

    if (ax <= 0) {
        return std::numeric_limits<int32_t>::min();
    }

    float x = static_cast<float>(ax) / static_cast<float>(1u << frac_bits);
    float result = 1.0f / std::sqrt(x);
    float scale = static_cast<float>(1u << frac_bits);
    long long result_scaled = llroundf(result * scale);

    return sat_cast<int32_t>(result_scaled);
}

// Forward to Priority 1 when NOT 32-bit
template<int Xb,
         typename std::enable_if<!is_rsqrt_32<Xb>::value, int>::type = 0>
inline typename StorageForBits<Xb>::type
xtensa_rsqrt_impl(typename StorageForBits<Xb>::type ax,
                  int frac_bits,
                  priority_tag<2>)
{
    return xtensa_rsqrt_impl<Xb>(ax, frac_bits, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
