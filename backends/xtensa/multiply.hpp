#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Multiply Implementation with Priority Dispatch
// ============================================================================
//
// This implementation uses priority_tag dispatch to select specialized
// multiplication paths based on operand/output bucket sizes:
//   - Priority 2: 8×8→8 (uses int32_t intermediate for NatureDSP vec_multiply calls)
//   - Priority 1: 16×16→16 (uses int32_t intermediate)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order (Priority 0 first,
// then Priority 1, then Priority 2) so that each level can call the next without forward declarations.

// -------- Priority 0: Generic Fallback → ReferenceBackend (define first) --------

template<int Xb, int Yb, int Ob>
inline typename StorageForBits<Ob>::type
xtensa_mul_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<0> )
{
    // Fallback to portable ReferenceBackend implementation
    return ReferenceBackend::template mul<Xb, Yb, Ob>(ax, by, out_frac_shift);
}

// -------- Priority 1: 16×16→16 Specialization --------

template<int Xb, int Yb, int Ob>
using is_16x16_to_16 = std::integral_constant<bool,
    IsBucket<Xb,16>::value && IsBucket<Yb,16>::value && IsBucket<Ob,16>::value>;

// Enabled when 16×16→16
template<int Xb, int Yb, int Ob,
         typename std::enable_if<is_16x16_to_16<Xb,Yb,Ob>::value, int>::type = 0>
inline int16_t
xtensa_mul_impl(int16_t ax, int16_t by, int out_frac_shift, priority_tag<1>)
{
    // TODO: Use NatureDSP function when available:
    // For now, use portable implementation with int32 intermediate
    int32_t prod = static_cast<int32_t>(ax) * static_cast<int32_t>(by);
    int32_t shifted = (out_frac_shift >= 0) ? round_shift(prod, out_frac_shift)
                                            : (prod << (-out_frac_shift));
    return sat_cast<int16_t>(shifted);
}

// Forward to Priority 0 when NOT 16×16→16
template<int Xb, int Yb, int Ob,
         typename std::enable_if<!is_16x16_to_16<Xb,Yb,Ob>::value, int>::type = 0>
inline typename StorageForBits<Ob>::type
xtensa_mul_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<1> )
{
    return xtensa_mul_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<0>{});
}

// -------- Priority 2: 8×8→8 Specialization --------

template<int Xb, int Yb, int Ob>
using is_8x8_to_8 = std::integral_constant<bool,
    IsBucket<Xb,8>::value && IsBucket<Yb,8>::value && IsBucket<Ob,8>::value>;

// Enabled when 8×8→8
template<int Xb, int Yb, int Ob,
         typename std::enable_if<is_8x8_to_8<Xb,Yb,Ob>::value, int>::type = 0>
inline int8_t
xtensa_mul_impl(int8_t ax, int8_t by, int out_frac_shift, priority_tag<2>)
{
    // TODO: Use NatureDSP function when available:
    // For now, use portable implementation with int32 intermediate
    int32_t prod = static_cast<int32_t>(ax) * static_cast<int32_t>(by);
    int32_t shifted = (out_frac_shift >= 0) ? round_shift(prod, out_frac_shift)
                                            : (prod << (-out_frac_shift));
    return sat_cast<int8_t>(shifted);
}

// Forward to Priority 1 when NOT 8×8→8
template<int Xb, int Yb, int Ob,
         typename std::enable_if<!is_8x8_to_8<Xb,Yb,Ob>::value, int>::type = 0>
inline typename StorageForBits<Ob>::type
xtensa_mul_impl( typename StorageForBits<Xb>::type ax,
                 typename StorageForBits<Yb>::type by,
                 int out_frac_shift, priority_tag<2> )
{
    return xtensa_mul_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
