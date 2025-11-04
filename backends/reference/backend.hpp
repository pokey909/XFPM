#pragma once
#include "../../helpers.hpp"

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
 * ReferenceBackend: Portable C++ implementation of fixed-point operations
 *
 * This backend provides reference implementations using standard C++ and
 * helper functions from helpers.hpp. All operations are portable and do not
 * depend on any specific hardware or library.
 *
 * Each operation is implemented in a separate header file in the reference/
 * directory and included above.
 */
struct ReferenceBackend {
    // Multiply operation
    template<int Xb, int Yb, int Ob>
    static typename StorageForBits<Ob>::type
    mul( typename StorageForBits<Xb>::type ax,
         typename StorageForBits<Yb>::type by,
         int out_frac_shift )
    {
        return detail::reference_mul<Xb, Yb, Ob>(ax, by, out_frac_shift);
    }

    // Divide operation
    template<int Xb, int Yb, int Ob>
    static typename StorageForBits<Ob>::type
    div( typename StorageForBits<Xb>::type ax,
         typename StorageForBits<Yb>::type by,
         int out_frac_shift )
    {
        return detail::reference_div<Xb, Yb, Ob>(ax, by, out_frac_shift);
    }

    // Logarithm operations (input converted to Q16.15, output as Q6.25)
    template<int Xb>
    static int32_t log2(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::reference_log2<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static int32_t logn(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::reference_logn<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static int32_t log10(typename StorageForBits<Xb>::type ax, int frac_bits) {
        return detail::reference_log10<Xb>(ax, frac_bits);
    }

    // Antilogarithm operations (input as Q6.25, output as Q16.15)
    template<int Xb>
    static int32_t antilog2(Storage_t<Xb> ax, int frac_bits) {
        return detail::reference_antilog2<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static int32_t antilogn(Storage_t<Xb> ax, int frac_bits) {
        return detail::reference_antilogn<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static int32_t antilog10(Storage_t<Xb> ax, int frac_bits) {
        return detail::reference_antilog10<Xb>(ax, frac_bits);
    }

    // Power operation
    template<int Xb, int Yb>
    static typename StorageForBits<Xb>::type
    pow(typename StorageForBits<Xb>::type base,
        typename StorageForBits<Yb>::type exponent,
        int base_frac_bits,
        int exp_frac_bits)
    {
        return detail::reference_pow<Xb, Yb>(base, exponent, base_frac_bits, exp_frac_bits);
    }

    // Square root operation (returns same Q format as input)
    template<int Xb>
    static typename StorageForBits<Xb>::type
    sqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
    {
        return detail::reference_sqrt<Xb>(ax, frac_bits);
    }

    // Reciprocal square root operation (returns same Q format as input)
    template<int Xb>
    static typename StorageForBits<Xb>::type
    rsqrt(typename StorageForBits<Xb>::type ax, int frac_bits)
    {
        return detail::reference_rsqrt<Xb>(ax, frac_bits);
    }

    // Array Shift/Scale operations (in-place)
    template<int Xb>
    static void
    array_shift(Storage_t<Xb>* arr, size_t length, int shift_amount)
    {
        detail::reference_array_shift<Xb>(arr, length, shift_amount);
    }

    template<int Xb>
    static void
    array_scale(Storage_t<Xb>* arr, size_t length,
                Storage_t<Xb> scale_factor, int scale_frac_bits)
    {
        detail::reference_array_scale<Xb>(arr, length, scale_factor, scale_frac_bits);
    }

    // Array Min/Max operations
    template<int Xb>
    static Storage_t<Xb>
    array_min(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::reference_array_min<Xb>(arr, length);
    }

    template<int Xb>
    static Storage_t<Xb>
    array_max(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::reference_array_max<Xb>(arr, length);
    }

    // Vector operations
    template<int Xb>
    static Storage_t<Xb>
    dot_product(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits)
    {
        return detail::reference_dot_product<Xb>(arr1, arr2, length, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    array_sum(const Storage_t<Xb>* arr, size_t length)
    {
        return detail::reference_array_sum<Xb>(arr, length);
    }

    // Trigonometric operations
    template<int Xb>
    static Storage_t<Xb>
    sin(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_sin<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    cos(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_cos<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    tan(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_tan<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    atan(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_atan<Xb>(ax, frac_bits);
    }

    // Hyperbolic functions
    template<int Xb>
    static Storage_t<Xb>
    tanh(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_tanh<Xb>(ax, frac_bits);
    }

    // Activation functions
    template<int Xb>
    static Storage_t<Xb>
    sigmoid(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_sigmoid<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    relu(Storage_t<Xb> ax, int frac_bits)
    {
        return detail::reference_relu<Xb>(ax, frac_bits);
    }

    template<int Xb>
    static void
    softmax(const Storage_t<Xb>* input, Storage_t<Xb>* output,
            size_t length, int frac_bits)
    {
        detail::reference_softmax<Xb>(input, output, length, frac_bits);
    }

    // Element-wise vector operations
    template<int Xb>
    static void
    array_elemult(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                  Storage_t<Xb>* output, size_t length, int frac_bits)
    {
        detail::reference_array_elemult<Xb>(arr1, arr2, output, length, frac_bits);
    }

    template<int Xb>
    static void
    array_add(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
              Storage_t<Xb>* output, size_t length)
    {
        detail::reference_array_add<Xb>(arr1, arr2, output, length);
    }

    template<int Xb>
    static void
    array_sub(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
              Storage_t<Xb>* output, size_t length)
    {
        detail::reference_array_sub<Xb>(arr1, arr2, output, length);
    }

    // Statistical vector operations
    template<int Xb>
    static Storage_t<Xb>
    array_mean(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::reference_array_mean<Xb>(arr, length, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    array_rms(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::reference_array_rms<Xb>(arr, length, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    array_variance(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::reference_array_variance<Xb>(arr, length, frac_bits);
    }

    template<int Xb>
    static Storage_t<Xb>
    array_stddev(const Storage_t<Xb>* arr, size_t length, int frac_bits)
    {
        return detail::reference_array_stddev<Xb>(arr, length, frac_bits);
    }

    // Future operations will be added here as methods that forward to
    // implementations in their respective category header files
};

} // namespace fp
