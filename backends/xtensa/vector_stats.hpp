#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include <cstddef>
#include <cmath>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Vector Statistical Operations with Priority Dispatch
// ============================================================================
//
// NatureDSP functions available:
//   - vec_mean16x16, vec_mean32x32: Array mean
//   - vec_rms16x16, vec_rms32x32: Root mean square
//   - vec_var16x16, vec_var32x32: Variance
//   - vec_stddev16x16, vec_stddev32x32: Standard deviation
//   - vec_power16x16, vec_power32x32: Sum of squares
//   - Note: No 8-bit variants available
//
// Statistical operations use priority_tag dispatch:
//   - Priority 2: 32-bit and 8-bit arrays (NatureDSP optimized when available)
//   - Priority 1: 16-bit arrays (NatureDSP optimized when available)
//   - Priority 0: Generic fallback to ReferenceBackend

// ========== MEAN ==========

template<int Xb>
inline Storage_t<Xb>
xtensa_array_mean_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                       priority_tag<0>)
{
    return detail::reference_array_mean<Xb>(arr, length, frac_bits);
}

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

template<int Xb>
using is_8bit = std::integral_constant<bool, IsBucket<Xb, 8>::value>;

template<int Xb>
using is_32bit = std::integral_constant<bool, IsBucket<Xb, 32>::value>;

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_mean_impl(const int16_t* arr, size_t length, int frac_bits,
                       priority_tag<1>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_mean16x16
    // vec_mean16x16: int16_t vec_mean16x16(const int16_t *x, int N)
    // Returns mean of all elements with saturation
    return vec_mean16x16(arr, static_cast<int>(length));
}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_mean_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                       priority_tag<1>)
{
    return xtensa_array_mean_impl<Xb>(arr, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_array_mean_impl(const int32_t* arr, size_t length, int frac_bits,
                       priority_tag<2>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_mean32x32
    // vec_mean32x32: int32_t vec_mean32x32(const int32_t *x, int N)
    // Returns mean of all elements with saturation
    return vec_mean32x32(arr, static_cast<int>(length));
}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_mean_impl(const int8_t* arr, size_t length, int frac_bits,
                       priority_tag<2>)
{
    // NatureDSP doesn't have vec_mean8x8 - use reference implementation
    return xtensa_array_mean_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_mean_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                       priority_tag<2>)
{
    return xtensa_array_mean_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

// ========== RMS ==========

template<int Xb>
inline Storage_t<Xb>
xtensa_array_rms_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                      priority_tag<0>)
{
    return detail::reference_array_rms<Xb>(arr, length, frac_bits);
}

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_rms_impl(const int16_t* arr, size_t length, int frac_bits,
                      priority_tag<1>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_rms16x16
    // vec_rms16x16: int16_t vec_rms16x16(const int16_t *x, int N)
    // Returns RMS of all elements with saturation
    return vec_rms16x16(arr, static_cast<int>(length));
}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_rms_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                      priority_tag<1>)
{
    return xtensa_array_rms_impl<Xb>(arr, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_array_rms_impl(const int32_t* arr, size_t length, int frac_bits,
                      priority_tag<2>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_rms32x32
    // vec_rms32x32: int32_t vec_rms32x32(const int32_t *x, int N)
    // Returns RMS of all elements with saturation
    return vec_rms32x32(arr, static_cast<int>(length));
}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_rms_impl(const int8_t* arr, size_t length, int frac_bits,
                      priority_tag<2>)
{
    // NatureDSP doesn't have vec_rms8x8 - use reference implementation
    return xtensa_array_rms_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_rms_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                      priority_tag<2>)
{
    return xtensa_array_rms_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

// ========== VARIANCE ==========

template<int Xb>
inline Storage_t<Xb>
xtensa_array_variance_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                           priority_tag<0>)
{
    return detail::reference_array_variance<Xb>(arr, length, frac_bits);
}

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_variance_impl(const int16_t* arr, size_t length, int frac_bits,
                           priority_tag<1>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_var16x16
    // vec_var16x16: int16_t vec_var16x16(const int16_t *x, int N)
    // Returns variance of all elements with saturation
    return vec_var16x16(arr, static_cast<int>(length));

}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_variance_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                           priority_tag<1>)
{
    return xtensa_array_variance_impl<Xb>(arr, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_array_variance_impl(const int32_t* arr, size_t length, int frac_bits,
                           priority_tag<2>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_var32x32
    // vec_var32x32: int32_t vec_var32x32(const int32_t *x, int N)
    // Returns variance of all elements with saturation
    return vec_var32x32(arr, static_cast<int>(length));

}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_variance_impl(const int8_t* arr, size_t length, int frac_bits,
                           priority_tag<2>)
{
    // NatureDSP doesn't have vec_var8x8 - use reference implementation
    return xtensa_array_variance_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_variance_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                           priority_tag<2>)
{
    return xtensa_array_variance_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

// ========== STANDARD DEVIATION ==========

template<int Xb>
inline Storage_t<Xb>
xtensa_array_stddev_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                         priority_tag<0>)
{
    return detail::reference_array_stddev<Xb>(arr, length, frac_bits);
}

template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_stddev_impl(const int16_t* arr, size_t length, int frac_bits,
                         priority_tag<1>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_stddev16x16
    // vec_stddev16x16: int16_t vec_stddev16x16(const int16_t *x, int N)
    // Returns standard deviation of all elements with saturation
    return vec_stddev16x16(arr, static_cast<int>(length));

}

template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_stddev_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                         priority_tag<1>)
{
    return xtensa_array_stddev_impl<Xb>(arr, length, frac_bits, priority_tag<0>{});
}

template<int Xb, EnableIf<is_32bit<Xb>::value> = 0>
inline int32_t
xtensa_array_stddev_impl(const int32_t* arr, size_t length, int frac_bits,
                         priority_tag<2>)
{
    if (length == 0) return 0;

    // Use NatureDSP vec_stddev32x32
    // vec_stddev32x32: int32_t vec_stddev32x32(const int32_t *x, int N)
    // Returns standard deviation of all elements with saturation
    return vec_stddev32x32(arr, static_cast<int>(length));

}

template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_stddev_impl(const int8_t* arr, size_t length, int frac_bits,
                         priority_tag<2>)
{
    // NatureDSP doesn't have vec_stddev8x8 - use reference implementation
    return xtensa_array_stddev_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

template<int Xb, EnableIf<!is_32bit<Xb>::value && !is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_stddev_impl(const Storage_t<Xb>* arr, size_t length, int frac_bits,
                         priority_tag<2>)
{
    return xtensa_array_stddev_impl<Xb>(arr, length, frac_bits, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
