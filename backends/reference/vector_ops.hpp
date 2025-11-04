#pragma once
#include "../../helpers.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Reference Vector Operations Implementation
// ============================================================================
//
// Array operations for dot product and sum.
// These operations work on arrays of fixed-point values stored as integers.

// Compute dot product of two arrays
// Performs fixed-point multiplication with proper scaling
template<int Xb>
inline Storage_t<Xb>
reference_dot_product(const Storage_t<Xb>* arr1, const Storage_t<Xb>* arr2, size_t length, int frac_bits)
{
    if (length == 0) {
        return 0;  // Return 0 for empty array
    }

    long long result = 0;
    for (size_t i = 0; i < length; ++i) {
        // Fixed-point multiply: multiply then shift right by frac_bits
        long long product = static_cast<long long>(arr1[i]) * static_cast<long long>(arr2[i]);
        result += round_shift(product, frac_bits);
    }
    return sat_cast<Storage_t<Xb>>(result);
}

// Compute sum of all elements in array
template<int Xb>
inline Storage_t<Xb>
reference_array_sum(const Storage_t<Xb>* arr, size_t length)
{
    if (length == 0) {
        return 0;  // Return 0 for empty array
    }

    Storage_t<Xb> sum = 0;
    for (size_t i = 0; i < length; ++i) {
        sum += arr[i];
    }
    return sum;
}

} // namespace detail
} // namespace fp
