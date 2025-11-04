#pragma once
#include "../../helpers.hpp"

namespace fp {

// Forward declare the backend struct
struct ReferenceBackend;

// Multiply operation implementation for ReferenceBackend
namespace detail {

template<int Xb, int Yb, int Ob>
inline typename StorageForBits<Ob>::type
reference_mul( typename StorageForBits<Xb>::type ax,
               typename StorageForBits<Yb>::type by,
               int out_frac_shift )
{
    using Out = typename StorageForBits<Ob>::type;
    long long prod = static_cast<long long>(ax) * static_cast<long long>(by);
    long long shifted = (out_frac_shift >= 0) ? round_shift(prod, out_frac_shift)
                                              : (prod << (-out_frac_shift));
    return sat_cast<Out>(shifted);
}

} // namespace detail
} // namespace fp
