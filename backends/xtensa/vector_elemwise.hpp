#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include "vector_helpers.hpp"
#include <cstddef>
#include <cassert>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Element-wise Vector Operations with Priority Dispatch
// ============================================================================
//
// Element-wise operations use priority_tag dispatch:
//   - Priority 2: 32-bit and 8-bit arrays (NatureDSP optimized when available)
//   - Priority 1: 16-bit arrays (NatureDSP optimized when available)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// NatureDSP functions available:
//   - vec_elemult16x16, vec_elemult32x32: Element-wise multiplication (combined with vec_shift for Q-format)
//   - vec_add16x16, vec_add32x32: Element-wise addition
//   - vec_add16x16_fast, vec_add32x32_fast: Optimized addition (requires 8-byte alignment, N%4==0)
//   - vec_elesub16x16, vec_elesub32x32: Element-wise subtraction
//   - vec_shift16x16, vec_shift32x32: Used for Q-format correction after elemult
//   - Note: No 8-bit variants available

// ========== ELEMENT-WISE MULTIPLY ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_array_elemult_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                          Storage_t<Xb>* output, size_t length, int frac_bits,
                          priority_tag<0>)
{
    return detail::reference_array_elemult<Xb>(arr1, arr2, output, length, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_array_elemult_impl(const int16_t* arr1, const int16_t* arr2,
                          int16_t* output, size_t length, int frac_bits,
                          priority_tag<1>)
{
    // NatureDSP vec_elemult16x16 performs simple 16x16 -> 16-bit multiplication with saturation.
    // For Q1.15 inputs, the intermediate product (e.g., 0.5*0.5 = 16384*16384 = 268M) vastly
    // exceeds 16-bit range and saturates. This is incompatible with Q-format arithmetic which
    // requires the full 32-bit product to extract the correct fractional bits.
    // Fall back to reference implementation which does proper Q-format multiply with rounding.
    return xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_array_elemult_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                          Storage_t<Xb>* output, size_t length, int frac_bits,
                          priority_tag<1>)
{
    return xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

template<int Xb>
using is_8bit = std::integral_constant<bool, IsBucket<Xb, 8>::value>;

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline void
xtensa_array_elemult_impl(const int32_t* arr1, const int32_t* arr2,
                          int32_t* output, size_t length, int frac_bits,
                          priority_tag<2>)
{
    // NatureDSP vec_elemult32x32 has the same issue as vec_elemult16x16:
    // it performs simple 32x32 -> 32-bit multiplication with saturation, which causes
    // the intermediate 64-bit product to overflow and saturate. Q-format arithmetic
    // requires the full 64-bit product to extract the correct fractional bits.
    // Fall back to reference implementation which does proper Q-format multiply with rounding.
    return xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline void
xtensa_array_elemult_impl(const int8_t* arr1, const int8_t* arr2,
                          int8_t* output, size_t length, int frac_bits,
                          priority_tag<2>)
{
    // NatureDSP doesn't have vec_elemult8x8 - use reference implementation
    return xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline void
xtensa_array_elemult_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                          Storage_t<Xb>* output, size_t length, int frac_bits,
                          priority_tag<2>)
{
    return xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<1>{});
}

// ========== ELEMENT-WISE ADD ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_array_add_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<0>)
{
    return detail::reference_array_add<Xb>(arr1, arr2, output, length);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_array_add_impl(const int16_t* arr1, const int16_t* arr2,
                      int16_t* output, size_t length, priority_tag<1>)
{
    // Use NatureDSP vec_add16x16
    // vec_add16x16: void vec_add16x16(int16_t *z, const int16_t *x, const int16_t *y, int N)
    // vec_add16x16_fast: void vec_add16x16_fast(int16_t *z, const int16_t *x, const int16_t *y, int N)
    //   - Requires 8-byte alignment and N%4==0
    // Performs: z[i] = x[i] + y[i] with saturation
    if (can_use_fast_variant(arr1, arr2, output, length)) {
        // Use fast variant
        vec_add16x16_fast(output, arr1, arr2, static_cast<int>(length));
    } else {
        // Use normal variant
        vec_add16x16(output, arr1, arr2, static_cast<int>(length));
    }
}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_array_add_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<1>)
{
    return xtensa_array_add_impl<Xb>(arr1, arr2, output, length, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline void
xtensa_array_add_impl(const int32_t* arr1, const int32_t* arr2,
                      int32_t* output, size_t length, priority_tag<2>)
{
    // Use NatureDSP vec_add32x32
    // vec_add32x32: void vec_add32x32(int32_t *z, const int32_t *x, const int32_t *y, int N)
    // vec_add32x32_fast: void vec_add32x32_fast(int32_t *z, const int32_t *x, const int32_t *y, int N)
    //   - Requires 8-byte alignment and N%4==0
    // Performs: z[i] = x[i] + y[i] with saturation
    if (can_use_fast_variant(arr1, arr2, output, length)) {
        // Use fast variant
        vec_add32x32_fast(output, arr1, arr2, static_cast<int>(length));
    } else {
        // Use normal variant
        vec_add32x32(output, arr1, arr2, static_cast<int>(length));
    }
}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline void
xtensa_array_add_impl(const int8_t* arr1, const int8_t* arr2,
                      int8_t* output, size_t length, priority_tag<2>)
{
    // NatureDSP doesn't have vec_add8x8 - use reference implementation
    return xtensa_array_add_impl<Xb>(arr1, arr2, output, length, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline void
xtensa_array_add_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<2>)
{
    return xtensa_array_add_impl<Xb>(arr1, arr2, output, length, priority_tag<1>{});
}

// ========== ELEMENT-WISE SUBTRACT ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_array_sub_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<0>)
{
    return detail::reference_array_sub<Xb>(arr1, arr2, output, length);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_array_sub_impl(const int16_t* arr1, const int16_t* arr2,
                      int16_t* output, size_t length, priority_tag<1>)
{
    // Use NatureDSP vec_elesub16x16
    // vec_elesub16x16: void vec_elesub16x16(int16_t *z, int16_t *x, int16_t *y, int N)
    // Performs: z[i] = x[i] - y[i] with saturation
    vec_elesub16x16(output, const_cast<int16_t*>(arr1), const_cast<int16_t*>(arr2), static_cast<int>(length));
}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_array_sub_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<1>)
{
    return xtensa_array_sub_impl<Xb>(arr1, arr2, output, length, priority_tag<0>{});
}

// -------- Priority 2: 32-bit and 8-bit Specializations --------

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline void
xtensa_array_sub_impl(const int32_t* arr1, const int32_t* arr2,
                      int32_t* output, size_t length, priority_tag<2>)
{
    // Use NatureDSP vec_elesub32x32
    // vec_elesub32x32: void vec_elesub32x32(int32_t *z, int32_t *x, int32_t *y, int N)
    // Performs: z[i] = x[i] - y[i] with saturation
    vec_elesub32x32(output, const_cast<int32_t*>(arr1), const_cast<int32_t*>(arr2), static_cast<int>(length));
}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline void
xtensa_array_sub_impl(const int8_t* arr1, const int8_t* arr2,
                      int8_t* output, size_t length, priority_tag<2>)
{
    // NatureDSP doesn't have vec_elesub8x8 - use reference implementation
    return xtensa_array_sub_impl<Xb>(arr1, arr2, output, length, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline void
xtensa_array_sub_impl(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                      Storage_t<Xb>* output, size_t length, priority_tag<2>)
{
    return xtensa_array_sub_impl<Xb>(arr1, arr2, output, length, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
