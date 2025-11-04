#pragma once
#include "../../helpers.hpp"
#include "../reference/backend.hpp"
#include <cstddef>

namespace fp {
namespace detail {

// ============================================================================
// Xtensa Vector Min/Max Implementation with Priority Dispatch
// ============================================================================
//
// Array min/max operations use priority_tag dispatch:
//   - Priority 2: 8-bit arrays (could use NatureDSP vec_min8/max8)
//   - Priority 1: 16-bit arrays (could use NatureDSP vec_min16/max16)
//   - Priority 0: Generic fallback to ReferenceBackend
//
// Note: All implementations must be defined in reverse priority order.

// ========== ARRAY MIN ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_array_min_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<0>)
{
    return ReferenceBackend::template array_min<Xb>(arr, length);
}

// -------- Priority 1: 16-bit Specialization --------

template<int Xb>
using is_16bit = std::integral_constant<bool, IsBucket<Xb, 16>::value>;

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_min_impl(const int16_t* arr, size_t length, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_min16 when available
    // For now, use portable implementation
    return xtensa_array_min_impl<Xb>(arr, length, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_min_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<1>)
{
    return xtensa_array_min_impl<Xb>(arr, length, priority_tag<0>{});
}

// -------- Priority 2: 8-bit Specialization --------

template<int Xb>
using is_8bit = std::integral_constant<bool, IsBucket<Xb, 8>::value>;

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_min_impl(const int8_t* arr, size_t length, priority_tag<2>)
{
    // TODO: Use NatureDSP vec_min8 when available
    // For now, use portable implementation
    return xtensa_array_min_impl<Xb>(arr, length, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 8-bit
template<int Xb, EnableIf<!is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_min_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<2>)
{
    return xtensa_array_min_impl<Xb>(arr, length, priority_tag<1>{});
}

// ========== ARRAY MAX ==========

// -------- Priority 0: Generic Fallback → ReferenceBackend --------

template<int Xb>
inline Storage_t<Xb>
xtensa_array_max_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<0>)
{
    return ReferenceBackend::template array_max<Xb>(arr, length);
}

// -------- Priority 1: 16-bit Specialization --------

// Enabled when 16-bit
template<int Xb, EnableIf<is_16bit<Xb>::value> = 0>
inline int16_t
xtensa_array_max_impl(const int16_t* arr, size_t length, priority_tag<1>)
{
    // TODO: Use NatureDSP vec_max16 when available
    // For now, use portable implementation
    return xtensa_array_max_impl<Xb>(arr, length, priority_tag<0>{});
}

// Forward to Priority 0 when NOT 16-bit
template<int Xb, EnableIf<!is_16bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_max_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<1>)
{
    return xtensa_array_max_impl<Xb>(arr, length, priority_tag<0>{});
}

// -------- Priority 2: 8-bit Specialization --------

// Enabled when 8-bit
template<int Xb, EnableIf<is_8bit<Xb>::value> = 0>
inline int8_t
xtensa_array_max_impl(const int8_t* arr, size_t length, priority_tag<2>)
{
    // TODO: Use NatureDSP vec_max8 when available
    // For now, use portable implementation
    return xtensa_array_max_impl<Xb>(arr, length, priority_tag<1>{});
}

// Forward to Priority 1 when NOT 8-bit
template<int Xb, EnableIf<!is_8bit<Xb>::value> = 0>
inline Storage_t<Xb>
xtensa_array_max_impl(const Storage_t<Xb>* arr, size_t length, priority_tag<2>)
{
    return xtensa_array_max_impl<Xb>(arr, length, priority_tag<1>{});
}

} // namespace detail
} // namespace fp
