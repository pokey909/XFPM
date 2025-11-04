#pragma once
#include <cstddef>
#include <cstdint>

namespace fp {
namespace detail {

// ============================================================================
// Helper Functions for NatureDSP Fast Variants
// ============================================================================

// Check if pointer is 8-byte aligned and length is multiple of 4 (for *_fast variants)
// Two-pointer version (for operations like dot_product)
template<typename T>
inline bool can_use_fast_variant(const T* ptr1, const T* ptr2, size_t length) {
    if (length < 4 || (length % 4) != 0) return false;
    return ((reinterpret_cast<uintptr_t>(ptr1) & 7) == 0) &&
           ((reinterpret_cast<uintptr_t>(ptr2) & 7) == 0);
}

// Three-pointer version (for operations like add with output)
template<typename T>
inline bool can_use_fast_variant(const T* ptr1, const T* ptr2, T* ptr3, size_t length) {
    if (length < 4 || (length % 4) != 0) return false;
    return ((reinterpret_cast<uintptr_t>(ptr1) & 7) == 0) &&
           ((reinterpret_cast<uintptr_t>(ptr2) & 7) == 0) &&
           ((reinterpret_cast<uintptr_t>(ptr3) & 7) == 0);
}

// Single-pointer version (for in-place operations like shift/scale)
template<typename T>
inline bool can_use_fast_variant(T* ptr, size_t length) {
    if (length < 4 || (length % 4) != 0) return false;
    return ((reinterpret_cast<uintptr_t>(ptr) & 7) == 0);
}

} // namespace detail
} // namespace fp
