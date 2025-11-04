#pragma once
#include "../../helpers.hpp"
#include <cmath>
#include <algorithm>
#include <cstddef>

namespace fp {
struct ReferenceBackend;

namespace detail {

// ============================================================================
// Reference Activation Functions Implementation
// ============================================================================
//
// Standard activation functions for neural networks:
//   - sigmoid: 1/(1 + e^-x), output range [0,1]
//   - relu: max(0, x)
//   - softmax: exp(x_i) / sum(exp(x_j)), array operation
//
// These functions work with fixed-point values and use the same Q format
// for input and output to maintain consistency with the storage type.

// ========== SIGMOID ==========
// Sigmoid: 1/(1 + e^-x)
// Input and output use same Q format
template<int Xb>
inline Storage_t<Xb>
reference_sigmoid(Storage_t<Xb> ax, int frac_bits)
{
    // Convert fixed-point to float
    float scale = static_cast<float>(1u << frac_bits);
    float x = static_cast<float>(ax) / scale;

    // Compute sigmoid: 1/(1 + e^-x)
    float result = 1.0f / (1.0f + std::exp(-x));

    // Convert back to fixed-point with same Q format
    long long result_scaled = llroundf(result * scale);

    return sat_cast<Storage_t<Xb>>(result_scaled);
}

// ========== RELU ==========
// ReLU: max(0, x)
// Input and output use same Q format
template<int Xb>
inline Storage_t<Xb>
reference_relu(Storage_t<Xb> ax, int frac_bits)
{
    // ReLU is simply max(0, x) in fixed-point
    return std::max(Storage_t<Xb>(0), ax);
}

// ========== SOFTMAX ==========
// Softmax: exp(x_i) / sum(exp(x_j))
// Array operation - takes input array, writes to output array
// Input and output use same Q format
template<int Xb>
inline void
reference_softmax(const Storage_t<Xb>* input, Storage_t<Xb>* output, size_t length, int frac_bits)
{
    if (length == 0) {
        return;  // Nothing to do for empty array
    }

    float scale = static_cast<float>(1u << frac_bits);

    // For numerical stability, subtract max value from all inputs
    Storage_t<Xb> max_val = input[0];
    for (size_t i = 1; i < length; ++i) {
        max_val = std::max(max_val, input[i]);
    }

    // Compute exp(x_i - max) for all elements and accumulate sum
    float sum = 0.0f;
    float* exp_vals = new float[length];

    for (size_t i = 0; i < length; ++i) {
        float x = static_cast<float>(input[i] - max_val) / scale;
        exp_vals[i] = std::exp(x);
        sum += exp_vals[i];
    }

    // Normalize by sum and convert back to fixed-point
    for (size_t i = 0; i < length; ++i) {
        float result = exp_vals[i] / sum;
        long long result_scaled = llroundf(result * scale);
        output[i] = sat_cast<Storage_t<Xb>>(result_scaled);
    }

    delete[] exp_vals;
}

} // namespace detail
} // namespace fp
