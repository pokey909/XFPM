#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include "vector_helpers.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Vector Scale Implementation with Priority Dispatch
// ============================================================================
//
// NatureDSP functions used:
//   - vec_shift16x16, vec_shift32x32: Arithmetic bit shift with saturation
//   - vec_shift16x16_fast, vec_shift32x32_fast: Fast variants (8-byte aligned, N%4==0)
//   - vec_scale16x16, vec_scale32x32: Scalar multiplication with saturation
//   - vec_scale16x16_fast, vec_scale32x32_fast: Fast variants (8-byte aligned, N%4==0)
//   - Note: No 8-bit variants available
//
// Array scale operations use priority_tag dispatch:
//   - Priority 2: 32-bit and 8-bit arrays (NatureDSP optimized when available)
//   - Priority 1: 16-bit arrays (NatureDSP optimized when available)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== ARRAY SHIFT ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_array_shift_impl(Storage_t<Xb>* arr, size_t length, int shift_amount, priority_tag<0>)
{
    return detail::reference_array_shift<Xb>(arr, length, shift_amount);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

template<int Xb>
using is_8bit = std::integral_constant<bool, IsBucket<Xb, 8>::value>;

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_array_shift_impl(int16_t* arr, size_t length, int shift_amount, priority_tag<1>)
{
    if (shift_amount == 0) return;  // No shift needed

    // NatureDSP vec_shift16x16: void vec_shift16x16(int16_t *y, const int16_t *x, int t, int N)
    // In-place operation: use same pointer for input and output
    // vec_shift16x16_fast: requires 8-byte alignment and N % 4 == 0
    if (can_use_fast_variant(arr, length)) {
        vec_shift16x16_fast(arr, arr, shift_amount, static_cast<int>(length));
    } else {
        vec_shift16x16(arr, arr, shift_amount, static_cast<int>(length));
    }
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_array_shift_impl(Storage_t<Xb>* arr, size_t length, int shift_amount, priority_tag<1>)
{
    return xtensa_array_shift_impl<Xb>(arr, length, shift_amount, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline void
xtensa_array_shift_impl(int32_t* arr, size_t length, int shift_amount, priority_tag<2>)
{
    if (shift_amount == 0) return;  // No shift needed

    // NatureDSP vec_shift32x32: void vec_shift32x32(int32_t *y, const int32_t *x, int t, int N)
    // In-place operation: use same pointer for input and output
    // vec_shift32x32_fast: requires 8-byte alignment and N % 4 == 0
    if (can_use_fast_variant(arr, length)) {
        vec_shift32x32_fast(arr, arr, shift_amount, static_cast<int>(length));
    } else {
        vec_shift32x32(arr, arr, shift_amount, static_cast<int>(length));
    }
}

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline void
xtensa_array_shift_impl(int8_t* arr, size_t length, int shift_amount, priority_tag<2>)
{
    // NatureDSP doesn't have vec_shift8x8 - use reference implementation
    return xtensa_array_shift_impl<Xb>(arr, length, shift_amount, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 32-bit and NOT 8-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline void
xtensa_array_shift_impl(Storage_t<Xb>* arr, size_t length, int shift_amount, priority_tag<2>)
{
    return xtensa_array_shift_impl<Xb>(arr, length, shift_amount, priority_tag<1>{});
}

// ========== ARRAY SCALE ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_array_scale_impl(Storage_t<Xb>* arr, size_t length,
                        Storage_t<Xb> scale_factor, int scale_frac_bits,
                        priority_tag<0>)
{
    return detail::reference_array_scale<Xb>(arr, length, scale_factor, scale_frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_array_scale_impl(int16_t* arr, size_t length,
                        int16_t scale_factor, int scale_frac_bits,
                        priority_tag<1>)
{
    // NatureDSP vec_scale16x16: void vec_scale16x16(int16_t *y, const int16_t *x, int16_t s, int N)
    // Performs: y[i] = sat(x[i] * s)
    // Then we need to shift right by scale_frac_bits to correct Q-format

    // Step 1: Multiply by scalar (in-place)
    if (can_use_fast_variant(arr, length)) {
        vec_scale16x16_fast(arr, arr, scale_factor, static_cast<int>(length));
    } else {
        vec_scale16x16(arr, arr, scale_factor, static_cast<int>(length));
    }

    // Step 2: Shift right by scale_frac_bits to correct Q-format (in-place)
    if (scale_frac_bits != 0) {
        if (can_use_fast_variant(arr, length)) {
            vec_shift16x16_fast(arr, arr, -scale_frac_bits, static_cast<int>(length));
        } else {
            vec_shift16x16(arr, arr, -scale_frac_bits, static_cast<int>(length));
        }
    }
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_array_scale_impl(Storage_t<Xb>* arr, size_t length,
                        Storage_t<Xb> scale_factor, int scale_frac_bits,
                        priority_tag<1>)
{
    return xtensa_array_scale_impl<Xb>(arr, length, scale_factor, scale_frac_bits, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

// Enabled when 32-bit
template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline void
xtensa_array_scale_impl(int32_t* arr, size_t length,
                        int32_t scale_factor, int scale_frac_bits,
                        priority_tag<2>)
{
    // NatureDSP vec_scale32x32: void vec_scale32x32(int32_t *y, const int32_t *x, int32_t s, int N)
    // Performs: y[i] = sat(x[i] * s)
    // Then we need to shift right by scale_frac_bits to correct Q-format

    // Step 1: Multiply by scalar (in-place)
    if (can_use_fast_variant(arr, length)) {
        vec_scale32x32_fast(arr, arr, scale_factor, static_cast<int>(length));
    } else {
        vec_scale32x32(arr, arr, scale_factor, static_cast<int>(length));
    }

    // Step 2: Shift right by scale_frac_bits to correct Q-format (in-place)
    if (scale_frac_bits != 0) {
        if (can_use_fast_variant(arr, length)) {
            vec_shift32x32_fast(arr, arr, -scale_frac_bits, static_cast<int>(length));
        } else {
            vec_shift32x32(arr, arr, -scale_frac_bits, static_cast<int>(length));
        }
    }
}

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline void
xtensa_array_scale_impl(int8_t* arr, size_t length,
                        int8_t scale_factor, int scale_frac_bits,
                        priority_tag<2>)
{
    // NatureDSP doesn't have vec_scale8x8 - use reference implementation
    return xtensa_array_scale_impl<Xb>(arr, length, scale_factor, scale_frac_bits, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 32-bit and NOT 8-bit
template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline void
xtensa_array_scale_impl(Storage_t<Xb>* arr, size_t length,
                        Storage_t<Xb> scale_factor, int scale_frac_bits,
                        priority_tag<2>)
{
    return xtensa_array_scale_impl<Xb>(arr, length, scale_factor, scale_frac_bits, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
