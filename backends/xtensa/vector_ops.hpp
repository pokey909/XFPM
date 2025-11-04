#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include "vector_helpers.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Vector Operations Implementation with Priority Dispatch
// ============================================================================
//
// Array operations use priority_tag dispatch:
//   - Priority 2: 32-bit and 8-bit arrays (NatureDSP optimized when available)
//   - Priority 1: 16-bit arrays (NatureDSP optimized when available)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// NatureDSP functions used:
//   - vec_dot16x16, vec_dot32x32: Dot product (returns int64_t accumulator)
//   - vec_dot16x16_fast, vec_dot32x32_fast: Optimized dot product (requires 8-byte alignment, N%4==0)
//   - vec_sum16x16, vec_sum32x32: Array sum (returns int64_t accumulator)
//   - vec_add16x16, vec_add32x32: Element-wise addition
//   - vec_power16x16, vec_power32x32: Sum of squares (for RMS)
//   - Note: No 8-bit variants available (vec_dot8x8, vec_sum8x8, etc.)
//   - Note: vec_dot16x16_fast returns int32_t (not int64_t)
//
// Note: All implementations must be defined in reverse priority order.

// ========== DOT PRODUCT ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_dot_product_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template dot_product<Xb>(arr1, arr2, length, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_dot_product_impl(const int16_t* arr1, const int16_t* arr2, size_t length, int frac_bits, priority_tag<1>)
{
    // Use NatureDSP vec_dot16x16 which returns 64-bit accumulator
    // vec_dot16x16: int64_t vec_dot16x16(const int16_t *x, const int16_t *y, int N)
    // vec_dot16x16_fast: int32_t vec_dot16x16_fast(const int16_t *x, const int16_t *y, int N)
    //   - Requires 8-byte alignment and N%4==0
    //   - Returns int32_t accumulator (may overflow for large vectors)
    // For QI.F format: Q(I,F) * Q(I,F) = Q(2I, 2F) in accumulator
    // We need to shift right by F to get back to Q(I,F)
    int64_t result64;
    if (can_use_fast_variant(arr1, arr2, length)) {
        // Use fast variant (returns int32_t, cast to int64_t)
        result64 = static_cast<int64_t>(vec_dot16x16_fast(arr1, arr2, static_cast<int>(length)));
    } else {
        // Use normal variant (returns int64_t)
        result64 = vec_dot16x16(arr1, arr2, static_cast<int>(length));
    }
    // Shift by frac_bits to convert from Q(2*F) back to Q(F)
    return sat_cast<int16_t>(round_shift(result64, frac_bits));
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_dot_product_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits, priority_tag<1>)
{
    return xtensa_dot_product_impl<Xb>(arr1, arr2, length, frac_bits, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

template<int Xb>
using is_8bit = std::integral_constant<bool, IsBucket<Xb, 8>::value>;

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_dot_product_impl(const int32_t* arr1, const int32_t* arr2, size_t length, int frac_bits, priority_tag<2>)
{
    // Use NatureDSP vec_dot32x32 which returns 64-bit accumulator
    // vec_dot32x32: int64_t vec_dot32x32(const int32_t *x, const int32_t *y, int N)
    // vec_dot32x32_fast: int64_t vec_dot32x32_fast(const int32_t *x, const int32_t *y, int N)
    //   - Requires 8-byte alignment and N%4==0
    // For QI.F format: Q(I,F) * Q(I,F) = Q(2I, 2F) in 64-bit accumulator
    // We need to shift right by F to get back to Q(I,F)
    int64_t result64;
    if (can_use_fast_variant(arr1, arr2, length)) {
        // Use fast variant
        result64 = vec_dot32x32_fast(arr1, arr2, static_cast<int>(length));
    } else {
        // Use normal variant
        result64 = vec_dot32x32(arr1, arr2, static_cast<int>(length));
    }
    // Shift by frac_bits to convert from Q(2*F) back to Q(F)
    return sat_cast<int32_t>(round_shift(result64, frac_bits));
}

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_dot_product_impl(const int8_t* arr1, const int8_t* arr2, size_t length, int frac_bits, priority_tag<2>)
{
    // NatureDSP doesn't have vec_dot8x8, use reference implementation
    return xtensa_dot_product_impl<Xb>(arr1, arr2, length, frac_bits, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 32-bit and NOT 8-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_dot_product_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits, priority_tag<2>)
{
    return xtensa_dot_product_impl<Xb>(arr1, arr2, length, frac_bits, priority_tag<1>{});
}

// ========== ARRAY SUM ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_array_sum_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<0>)
{
    return ReferenceBackend::template array_sum<Xb>(arr, length);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_sum_impl(const int16_t* arr, size_t length, priority_tag<1>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_sum16x16
    // vec_sum16x16: int64_t vec_sum16x16(const int16_t *x, int N)
    // Returns sum of all elements in 64-bit accumulator
    int64_t sum64 = vec_sum16x16(arr, static_cast<int>(length));
    return sat_cast<int16_t>(sum64);
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_sum_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<1>)
{
    return xtensa_array_sum_impl<Xb>(arr, length, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_array_sum_impl(const int32_t* arr, size_t length, priority_tag<2>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_sum32x32
    // vec_sum32x32: int64_t vec_sum32x32(const int32_t *x, int N)
    // Returns sum of all elements in 64-bit accumulator
    int64_t sum64 = vec_sum32x32(arr, static_cast<int>(length));
    return sat_cast<int32_t>(sum64);
}

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_sum_impl(const int8_t* arr, size_t length, priority_tag<2>)
{
    // NatureDSP doesn't have vec_sum8x8 - use reference implementation
    return xtensa_array_sum_impl<Xb>(arr, length, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 32-bit and NOT 8-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_sum_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<2>)
{
    return xtensa_array_sum_impl<Xb>(arr, length, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
