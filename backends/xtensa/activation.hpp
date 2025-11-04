#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Activation Functions Implementation with Priority Dispatch
// ============================================================================
//
// Activation functions use priority_tag dispatch:
//   - Priority 1: 16-bit operations (could use NatureDSP optimizations)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== SIGMOID ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_sigmoid_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template sigmoid<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_sigmoid_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP sigmoid optimization when available
    // For now, use portable implementation
    return xtensa_sigmoid_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_sigmoid_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_sigmoid_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== RELU ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_relu_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template relu<Xb>(ax, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_relu_impl(int16_t ax, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_max optimization when available
    // ReLU can be implemented as max(0, x) which NatureDSP supports efficiently
    // For now, use portable implementation
    return xtensa_relu_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_relu_impl(Storage_t<Xb> ax, int frac_bits, priority_tag<1>)
{
    return xtensa_relu_impl<Xb>(ax, frac_bits, priority_tag<0>{});
}

// ========== SOFTMAX ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline void
xtensa_softmax_impl(const Storage_t<Xb>* input, Storage_t<Xb>* output,
                    size_t length, int frac_bits, priority_tag<0>)
{
    return ReferenceBackend::template softmax<Xb>(input, output, length, frac_bits);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline void
xtensa_softmax_impl(const int16_t* input, int16_t* output,
                    size_t length, int frac_bits, priority_tag<1>)
{
    // TODO: Use NatureDSP vector operations for softmax when available
    // - vec_max16 for finding maximum
    // - vec_antilogn for computing exponentials
    // - vec_add16 for accumulation
    // For now, use portable implementation
    return xtensa_softmax_impl<Xb>(input, output, length, frac_bits, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline void
xtensa_softmax_impl(const Storage_t<Xb>* input, Storage_t<Xb>* output,
                    size_t length, int frac_bits, priority_tag<1>)
{
    return xtensa_softmax_impl<Xb>(input, output, length, frac_bits, priority_tag<0>{});
}

} // namespace detail
} // namespace fp
