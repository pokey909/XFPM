#pragma once
#include "../../helpers.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Reference Element-wise Vector Operations Implementation
// ============================================================================
//
// Element-wise operations on arrays of fixed-point values.
// These operations apply an operation to corresponding elements of arrays.

// Element-wise multiplication: output[i] = arr1[i] * arr2[i]
// Uses proper fixed-point multiply with shift
template<int Xb>
inline void
reference_array_elemult(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                        Storage_t<Xb>* output, size_t length, int frac_bits)
{
    if (length == 0) return;

    for (size_t i = 0; i < length; ++i) {
        // Use proper fixed-point multiply with rounding
        long long product = static_cast<long long>(arr1[i]) * static_cast<long long>(arr2[i]);
        output[i] = sat_cast<Storage_t<Xb>>(round_shift(product, frac_bits));
    }
}

// Element-wise addition: output[i] = arr1[i] + arr2[i]
template<int Xb>
inline void
reference_array_add(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                    Storage_t<Xb>* output, size_t length)
{
    if (length == 0) return;

    for (size_t i = 0; i < length; ++i) {
        long long sum = static_cast<long long>(arr1[i]) + static_cast<long long>(arr2[i]);
        output[i] = sat_cast<Storage_t<Xb>>(sum);
    }
}

// Element-wise subtraction: output[i] = arr1[i] - arr2[i]
template<int Xb>
inline void
reference_array_sub(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2,
                    Storage_t<Xb>* output, size_t length)
{
    if (length == 0) return;

    for (size_t i = 0; i < length; ++i) {
        long long diff = static_cast<long long>(arr1[i]) - static_cast<long long>(arr2[i]);
        output[i] = sat_cast<Storage_t<Xb>>(diff);
    }
}

} // namespace detail
} // namespace fp
