#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"

// Include NatureDSP library headers
// Note: These headers must be available in the include path when building for Xtensa
#ifdef __XTENSA__
#include <NatureDSP_Signal.h>
#include <NatureDSP_types.h>
#endif

// Include all category implementations
#include "multiply.hpp"
#include "divide.hpp"
#include "logarithm.hpp"
#include "antilogarithm.hpp"
#include "power.hpp"
#include "sqrt.hpp"
#include "rsqrt.hpp"
// Future category includes will go here:
#include "trigonometric.hpp"
#include "hyperbolic.hpp"
#include "activation.hpp"
#include "vector_ops.hpp"
#include "vector_scale.hpp"
#include "vector_minmax.hpp"
#include "vector_elemwise.hpp"
#include "vector_stats.hpp"

namespace fp {

/**
 * XtensaBackend: Optimized implementation using NatureDSP library
 *
 * This backend provides hardware-accelerated implementations of fixed-point
 * operations using the Cadence Xtensa HiFi4 Audio Engine NatureDSP library.
 *
 * The backend uses priority_tag dispatch to select the most specialized
 * implementation available:
 *  - Priority 2: 8×8→8 bit operations (int8_t operands/output)
 *  - Priority 1: 16×16→16 bit operations (int16_t operands/output)
 *  - Priority 0: Generic fallback to ReferenceBackend
 *
 * Each operation is implemented in a separate header file in the xtensa/
 * directory and uses SFINAE (enable_if) to select the appropriate path
 * based on the bucket sizes of the operands and output.
 *
 * When NatureDSP library is not available or for unsupported combinations,
 * operations fall back to the portable ReferenceBackend implementation.
 */
struct XtensaBackend {
    // Multiply operation with priority dispatch
    template<int Xb, int Yb, int Ob>
    static typename StorageForBits<Ob>::type
    mul( typename StorageForBits<Xb>::type ax,
         typename StorageForBits<Yb>::type by,
         int out_frac_shift )
    {
        return detail::xtensa_mul_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<2>{});
    }

    // Divide operation with priority dispatch
    template<int Xb, int Yb, int Ob>
    static typename StorageForBits<Ob>::type
    div( typename StorageForBits<Xb>::type ax,
         typename StorageForBits<Yb>::type by,
         int out_frac_shift )
    {
        return detail::xtensa_div_impl<Xb, Yb, Ob>(ax, by, out_frac_shift, priority_tag<2>{});
    }

    // Logarithm operations with priority dispatch
    template<int Xb>
    static int32_t log2(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::xtensa_log2_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static int32_t logn(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::xtensa_logn_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static int32_t log10(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::xtensa_log10_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    // Antilogarithm operations with priority dispatch
    template<int Xb>
    static int32_t antilog2(Storage_t<Xb> ax, int frac_bits) {
        return detail::xtensa_antilog2_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static int32_t antilogn(Storage_t<Xb> ax, int frac_bits) {
        return detail::xtensa_antilogn_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static int32_t antilog10(Storage_t<Xb> ax, int frac_bits) {
        return detail::xtensa_antilog10_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    // Power operation with priority dispatch
    template<int Xb, int Yb>
    static typename StorageForBits<Xb>::type
    pow(typename StorageForBits<Xb>::type base,
        typename StorageForBits<Yb>::type exponent,
        int base_frac_bits,
        int exp_frac_bits)
    {
        return detail::xtensa_pow_impl<Xb, Yb>(base, exponent, base_frac_bits, exp_frac_bits, priority_tag<1>{});
    }

    // Square root operation with priority dispatch
    template<int Xb>
    static typename StorageForBits<Xb>::type
    sqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
    {
        return detail::xtensa_sqrt_impl<Xb>(ax, frac_bits, priority_tag<2>{});
    }

    // Reciprocal square root operation with priority dispatch
    template<int Xb>
    static typename StorageForBits<Xb>::type
    rsqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
    {
        return detail::xtensa_rsqrt_impl<Xb>(ax, frac_bits, priority_tag<2>{});
    }

    // Array Shift/Scale operations with priority dispatch (in-place)
    template<int Xb>
    static void
    array_shift(Storage_t<Xb>* arr, size_t length, int shift_amount)
    {
        detail::xtensa_array_shift_impl<Xb>(arr, length, shift_amount, priority_tag<2>{});
    }

    template<int Xb>
    static void
    array_scale(Storage_t<Xb>* arr, size_t length,
                Storage_t<Xb> scale_factor, int scale_frac_bits)
    {
        detail::xtensa_array_scale_impl<Xb>(arr, length, scale_factor, scale_frac_bits, priority_tag<2>{});
    }

    // Array Min/Max operations with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    array_min(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::xtensa_array_min_impl<Xb>(arr, length, priority_tag<2>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    array_max(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::xtensa_array_max_impl<Xb>(arr, length, priority_tag<2>{});
    }

    // Vector operations with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    dot_product(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits)
    {
        return detail::xtensa_dot_product_impl<Xb>(arr1, arr2, length, frac_bits, priority_tag<2>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    array_sum(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::xtensa_array_sum_impl<Xb>(arr, length, priority_tag<2>{});
    }

    // Trigonometric operations with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    sin(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_sin_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    cos(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_cos_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    tan(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_tan_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    atan(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_atan_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    // Hyperbolic functions with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    tanh(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_tanh_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    // Activation functions with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    sigmoid(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_sigmoid_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    relu(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::xtensa_relu_impl<Xb>(ax, frac_bits, priority_tag<1>{});
    }

    template<int Xb>
    static void
    softmax(const Storage_t<Xb>* input, Storage_t<Xb>* output,
            size_t length, int frac_bits)
    {
        detail::xtensa_softmax_impl<Xb>(input, output, length, frac_bits, priority_tag<1>{});
    }

    // Element-wise vector operations with priority dispatch
    template<int Xb>
    static void
    array_elemult(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                  Storage_t<Xb>* output, size_t length, int frac_bits)
    {
        detail::xtensa_array_elemult_impl<Xb>(arr1, arr2, output, length, frac_bits, priority_tag<2>{});
    }

    template<int Xb>
    static void
    array_add(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
              Storage_t<Xb>* output, size_t length)
    {
        detail::xtensa_array_add_impl<Xb>(arr1, arr2, output, length, priority_tag<2>{});
    }

    template<int Xb>
    static void
    array_sub(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
              Storage_t<Xb>* output, size_t length)
    {
        detail::xtensa_array_sub_impl<Xb>(arr1, arr2, output, length, priority_tag<2>{});
    }

    // Statistical vector operations with priority dispatch
    template<int Xb>
    static Storage_t<Xb>
    array_mean(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::xtensa_array_mean_impl<Xb>(arr, length, frac_bits, priority_tag<2>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    array_rms(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::xtensa_array_rms_impl<Xb>(arr, length, frac_bits, priority_tag<2>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    array_variance(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::xtensa_array_variance_impl<Xb>(arr, length, frac_bits, priority_tag<2>{});
    }

    template<int Xb>
    static Storage_t<Xb>
    array_stddev(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::xtensa_array_stddev_impl<Xb>(arr, length, frac_bits, priority_tag<2>{});
    }

    // Future operations will be added here as methods that forward to
    // priority-dispatched implementations in their respective category header files
};

} // namespace fp
