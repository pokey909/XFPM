#pragma once
#include "../../helpers.hpp"
#include <algorithm>
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Reference Vector Min/Max Implementation
// ============================================================================
//
// Array operations that find minimum/maximum values in arrays.
// These operations work on arrays of fixed-point values stored as integers.

// Find minimum value in array
template<int Xb>
inline Storage_t<Xb>
reference_array_min(const Storage_t<Xb>* arr, size_t length)
{
    if (length == 0) {
        return 0;  // Return 0 for empty array
    }

    Storage_t<Xb> min_val = arr[0];
    for (size_t i = 1; i < length; ++i) {
        min_val = std::min(min_val, arr[i]);
    }
    return min_val;
}

// Find maximum value in array
template<int Xb>
inline Storage_t<Xb>
reference_array_max(const Storage_t<Xb>* arr, size_t length)
{
    if (length == 0) {
        return 0;  // Return 0 for empty array
    }

    Storage_t<Xb> max_val = arr[0];
    for (size_t i = 1; i < length; ++i) {
        max_val = std::max(max_val, arr[i]);
    }
    return max_val;
}

} // namespace detail
} // namespace fp
