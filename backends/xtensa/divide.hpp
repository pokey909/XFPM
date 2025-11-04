#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Divide Implementation with Priority Dispatch
// ============================================================================
//
// This implementation uses priority_tag dispatch to select specialized
// division paths based on operand/output bucket sizes:
//   - Priority 2: 8÷8→8 (uses int32_t intermediate for NatureDSP scl_divide16x16)
//   - Priority 1: 16÷16→16 (uses int64_t intermediate for NatureDSP scl_divide32x32)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order (Priority 0 first,
// then Priority 1, then Priority 2) so that each level can call the next without forward declarations.

// -------- Priority 0: Generic Fallback → ReferenceBackend (define first) --------

template<int Xb, int Yb, int Ob>
inline typename StorageForBits<Ob>::type
xtensa_div_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<0> )
{
    // Fallback to portable ReferenceBackend implementation
    return ReferenceBackend::template div<Xb, Yb, Ob>(ax, by, out_frac_shift);
}

// -------- Priority 1: 16÷16→16 Specialization --------

template<int Xb, int Yb, int Ob>
using is_16div16_to_16 = std::integral_constant<bool,
    IsBucket<Xb,16>::value && IsBucket<Yb,16>::value && IsBucket<Ob,16>::value>;

// Enabled when 16÷16→16
template<int Xb, int Yb, int Ob,
         typename std::enable_if<is_16div16_to_16<Xb,Yb,Ob>::value, int>::type = 0>
inline int16_t
xtensa_div_impl(int16_t ax, int16_t by, int out_frac_shift, priority_tag<1>)
{
    // TODO: Use NatureDSP scl_divide16x16 when available
    // For now, use portable implementation with int64 intermediate

    if (by == 0) {
        return (ax >= 0) ? std::numeric_limits<int16_t>::max()
                         : std::numeric_limits<int16_t>::min();
    }

    int shift = -out_frac_shift;
    int64_t dividend = static_cast<int64_t>(ax);
    int64_t divisor = static_cast<int64_t>(by);

    if (shift >= 0) {
        dividend = dividend << shift;
    } else {
        divisor = divisor << (-shift);
    }

    int64_t quotient;
    if ((dividend >= 0 && divisor > 0) || (dividend < 0 && divisor < 0)) {
        quotient = (dividend + divisor/2) / divisor;
    } else {
        quotient = (dividend - divisor/2) / divisor;
    }

    return sat_cast<int16_t>(quotient);
}

// Forward to Priority 0 when NOT 16÷16→16
template<int Xb, int Yb, int Ob,
         typename std::enable_if<!is_16div16_to_16<Xb,Yb,Ob>::value, int>::type = 0>
inline typename StorageForBits<Ob>::type
xtensa_div_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<1> )
{
    return xtensa_div_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<0>{});
}

// -------- Priority 2: 8÷8→8 Specialization --------

template<int Xb, int Yb, int Ob>
using is_8div8_to_8 = std::integral_constant<bool,
    IsBucket<Xb,8>::value && IsBucket<Yb,8>::value && IsBucket<Ob,8>::value>;

// Enabled when 8÷8→8
template<int Xb, int Yb, int Ob,
         typename std::enable_if<is_8div8_to_8<Xb,Yb,Ob>::value, int>::type = 0>
inline int8_t
xtensa_div_impl(int8_t ax, int8_t by, int out_frac_shift, priority_tag<2>)
{
    // TODO: Use NatureDSP scl_divide16x16 when available
    // For now, use portable implementation with int32 intermediate

    if (by == 0) {
        return (ax >= 0) ? std::numeric_limits<int8_t>::max()
                         : std::numeric_limits<int8_t>::min();
    }

    int shift = -out_frac_shift;
    int32_t dividend = static_cast<int32_t>(ax);
    int32_t divisor = static_cast<int32_t>(by);

    if (shift >= 0) {
        dividend = dividend << shift;
    } else {
        divisor = divisor << (-shift);
    }

    int32_t quotient;
    if ((dividend >= 0 && divisor > 0) || (dividend < 0 && divisor < 0)) {
        quotient = (dividend + divisor/2) / divisor;
    } else {
        quotient = (dividend - divisor/2) / divisor;
    }

    return sat_cast<int8_t>(quotient);
}

// Forward to Priority 1 when NOT 8÷8→8
template<int Xb, int Yb, int Ob,
         typename std::enable_if<!is_8div8_to_8<Xb,Yb,Ob>::value, int>::type = 0>
inline typename StorageForBits<Ob>::type
xtensa_div_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<2> )
{
    return xtensa_div_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
