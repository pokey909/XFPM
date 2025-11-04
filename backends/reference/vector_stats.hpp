#pragma once
#include "../../helpers.hpp"
#include <cstddef>
#include <cmath>

namespace fp {
namespace detail {

// ============================================================================
// Reference Vector Statistical Operations Implementation
// ============================================================================
//
// Statistical operations on arrays of fixed-point values.

// Mean: sum of elements / count
template<int Xb>
inline Storage_t<Xb>
reference_array_mean(const Storage_t<Xb>* arr, size_t length, int frac_bits)
{
    if (length == 0) return 0;

    // Sum all elements (using wider type to avoid overflow)
    long long sum = 0;
    for (size_t i = 0; i < length; ++i) {
        sum += arr[i];
    }

    // Divide by length (shift right by log2(length) for power-of-2, or actual divide)
    // For simplicity, use integer division
    long long mean = sum / static_cast<long long>(length);
    return sat_cast<Storage_t<Xb>>(mean);
}

// RMS (Root Mean Square): sqrt(sum(x^2) / N)
template<int Xb>
inline Storage_t<Xb>
reference_array_rms(const Storage_t<Xb>* arr, size_t length, int frac_bits)
{
    if (length == 0) return 0;

    // Sum of squares
    long long sum_squares = 0;
    for (size_t i = 0; i < length; ++i) {
        long long val = arr[i];
        long long square = (val * val) >> frac_bits;  // Fixed-point square with shift
        sum_squares += square;
    }

    // Mean of squares
    long long mean_square = sum_squares / static_cast<long long>(length);

    // Square root (convert to/from float for simplicity in reference implementation)
    float mean_sq_float = static_cast<float>(mean_square) / static_cast<float>(1u << frac_bits);
    float rms_float = std::sqrt(mean_sq_float);
    long long rms_scaled = llroundf(rms_float * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(rms_scaled);
}

// Variance: sum((x - mean)^2) / N
template<int Xb>
inline Storage_t<Xb>
reference_array_variance(const Storage_t<Xb>* arr, size_t length, int frac_bits)
{
    if (length == 0) return 0;

    // First compute mean
    Storage_t<Xb> mean = reference_array_mean<Xb>(arr, length, frac_bits);

    // Sum of squared deviations
    long long sum_sq_dev = 0;
    for (size_t i = 0; i < length; ++i) {
        long long deviation = static_cast<long long>(arr[i]) - static_cast<long long>(mean);
        long long sq_dev = (deviation * deviation) >> frac_bits;  // Fixed-point square
        sum_sq_dev += sq_dev;
    }

    // Divide by N
    long long variance = sum_sq_dev / static_cast<long long>(length);
    return sat_cast<Storage_t<Xb>>(variance);
}

// Standard deviation: sqrt(variance)
template<int Xb>
inline Storage_t<Xb>
reference_array_stddev(const Storage_t<Xb>* arr, size_t length, int frac_bits)
{
    if (length == 0) return 0;

    // Get variance
    Storage_t<Xb> variance = reference_array_variance<Xb>(arr, length, frac_bits);

    // Square root (convert to/from float for simplicity)
    float var_float = static_cast<float>(variance) / static_cast<float>(1u << frac_bits);
    float stddev_float = std::sqrt(var_float);
    long long stddev_scaled = llroundf(stddev_float * static_cast<float>(1u << frac_bits));

    return sat_cast<Storage_t<Xb>>(stddev_scaled);
}

} // namespace detail
} // namespace fp
