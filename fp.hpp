#pragma once
#include <type_traits>
#include <cstdint>
#include <cmath>
#include "helpers.hpp"     // must provide: StorageForBits<>, sat_cast<>
#include "backends/reference/backend.hpp"  // ReferenceBackend implementation
#ifdef __XTENSA__
#include "backends/xtensa/backend.hpp"     // XtensaBackend implementation (only for Xtensa builds)
#endif

namespace fp {

template<int I, int F, typename Backend = ReferenceBackend>
struct FixedPoint {
    static_assert(I >= 0 && F >= 0, "I and F must be non-negative");
    static constexpr int int_bits   = I;
    static constexpr int frac_bits  = F;
    static constexpr int total_bits = I + F;

    using storage_t    = typename StorageForBits< total_bits >::type;
    using backend_type = Backend;

    storage_t raw_;

    constexpr FixedPoint() : raw_(0) {}
    constexpr explicit FixedPoint(storage_t raw) : raw_(raw) {}

    // Round-to-nearest input quantization
    static FixedPoint from_float(float v) {
        const float scale = float(1u << F);
        long long q = llroundf(v * scale);
        return FixedPoint( fp::sat_cast<storage_t>(q) );
    }

    float to_float() const {
        return static_cast<float>(raw_) / static_cast<float>(1u << F);
    }

    storage_t raw() const { return raw_; }

    // Core compile-time routed multiply (explicit OUT_I/OUT_F)
    template<int OUT_I, int OUT_F, typename Other>
    auto mul(const Other& rhs) const {
        using Out = FixedPoint<OUT_I, OUT_F, Backend>;
        constexpr int Xb = total_bits;
        constexpr int Yb = Other::total_bits;
        constexpr int Ob = Out::total_bits;

        using Ax = typename StorageForBits<Xb>::type;
        using By = typename StorageForBits<Yb>::type;
        using Ro = typename StorageForBits<Ob>::type;

        constexpr int frac_in  = F + Other::frac_bits;
        constexpr int frac_out = OUT_F;
        constexpr int shift    = frac_in - frac_out;

        Ax ax = static_cast<Ax>(raw_);
        By by = static_cast<By>(rhs.raw());
        Ro ro = Backend::template mul<Xb,Yb,Ob>(ax, by, shift);
        return Out(ro);
    }

    // Ergonomic multiply: default to SAME Q as lhs (no ambiguity)
    template<typename Other>
    auto operator*(const Other& rhs) const {
        return this->template mul<I, F>(rhs);
    }

    // Core compile-time routed divide (explicit OUT_I/OUT_F)
    template<int OUT_I, int OUT_F, typename Other>
    auto div(const Other& rhs) const {
        using Out = FixedPoint<OUT_I, OUT_F, Backend>;
        constexpr int Xb = total_bits;
        constexpr int Yb = Other::total_bits;
        constexpr int Ob = Out::total_bits;

        using Ax = typename StorageForBits<Xb>::type;
        using By = typename StorageForBits<Yb>::type;
        using Ro = typename StorageForBits<Ob>::type;

        constexpr int frac_in  = F - Other::frac_bits;
        constexpr int frac_out = OUT_F;
        constexpr int shift    = frac_in - frac_out;

        Ax ax = static_cast<Ax>(raw_);
        By by = static_cast<By>(rhs.raw());
        Ro ro = Backend::template div<Xb,Yb,Ob>(ax, by, shift);
        return Out(ro);
    }

    // Ergonomic divide: default to SAME Q as lhs (no ambiguity)
    template<typename Other>
    auto operator/(const Other& rhs) const {
        return this->template div<I, F>(rhs);
    }

    // Core compile-time routed addition (explicit OUT_I/OUT_F)
    template<int OUT_I, int OUT_F, typename Other>
    auto add(const Other& rhs) const {
        using Out = FixedPoint<OUT_I, OUT_F, Backend>;
        constexpr int Xb = total_bits;
        constexpr int Yb = Other::total_bits;
        constexpr int Ob = Out::total_bits;

        using Ax = typename StorageForBits<Xb>::type;
        using By = typename StorageForBits<Yb>::type;
        using Ro = typename StorageForBits<Ob>::type;

        // Align operands to output fractional bits
        constexpr int shift_lhs = F - OUT_F;
        constexpr int shift_rhs = Other::frac_bits - OUT_F;

        Ax ax = static_cast<Ax>(raw_);
        By by = static_cast<By>(rhs.raw());

        // Shift operands to align fractional bits (with rounding)
        long long lhs_aligned = (shift_lhs == 0) ? ax : round_shift(ax, shift_lhs);
        long long rhs_aligned = (shift_rhs == 0) ? by : round_shift(by, shift_rhs);

        // Perform addition
        long long sum = lhs_aligned + rhs_aligned;

        // Saturate to output type
        Ro ro = fp::sat_cast<Ro>(sum);
        return Out(ro);
    }

    // Ergonomic addition: default to SAME Q as lhs
    template<typename Other>
    auto operator+(const Other& rhs) const {
        return this->template add<I, F>(rhs);
    }

    // Core compile-time routed subtraction (explicit OUT_I/OUT_F)
    template<int OUT_I, int OUT_F, typename Other>
    auto sub(const Other& rhs) const {
        using Out = FixedPoint<OUT_I, OUT_F, Backend>;
        constexpr int Xb = total_bits;
        constexpr int Yb = Other::total_bits;
        constexpr int Ob = Out::total_bits;

        using Ax = typename StorageForBits<Xb>::type;
        using By = typename StorageForBits<Yb>::type;
        using Ro = typename StorageForBits<Ob>::type;

        // Align operands to output fractional bits
        constexpr int shift_lhs = F - OUT_F;
        constexpr int shift_rhs = Other::frac_bits - OUT_F;

        Ax ax = static_cast<Ax>(raw_);
        By by = static_cast<By>(rhs.raw());

        // Shift operands to align fractional bits (with rounding)
        long long lhs_aligned = (shift_lhs == 0) ? ax : round_shift(ax, shift_lhs);
        long long rhs_aligned = (shift_rhs == 0) ? by : round_shift(by, shift_rhs);

        // Perform subtraction
        long long diff = lhs_aligned - rhs_aligned;

        // Saturate to output type
        Ro ro = fp::sat_cast<Ro>(diff);
        return Out(ro);
    }

    // Ergonomic subtraction: default to SAME Q as lhs
    template<typename Other>
    auto operator-(const Other& rhs) const {
        return this->template sub<I, F>(rhs);
    }

    // Comparison operators - handle mixed Q-format comparisons
    // by aligning to common fractional bits
    template<typename Other>
    bool operator<(const Other& rhs) const {
        // Align to max fractional bits for accurate comparison
        constexpr int max_frac = (F > Other::frac_bits) ? F : Other::frac_bits;
        constexpr int shift_lhs = F - max_frac;
        constexpr int shift_rhs = Other::frac_bits - max_frac;

        long long lhs_aligned = (shift_lhs == 0) ? raw_ : round_shift(raw_, shift_lhs);
        long long rhs_aligned = (shift_rhs == 0) ? rhs.raw() : round_shift(rhs.raw(), shift_rhs);

        return lhs_aligned < rhs_aligned;
    }

    template<typename Other>
    bool operator>(const Other& rhs) const {
        return rhs < *this;
    }

    template<typename Other>
    bool operator<=(const Other& rhs) const {
        return !(rhs < *this);
    }

    template<typename Other>
    bool operator>=(const Other& rhs) const {
        return !(*this < rhs);
    }

    template<typename Other>
    bool operator==(const Other& rhs) const {
        constexpr int max_frac = (F > Other::frac_bits) ? F : Other::frac_bits;
        constexpr int shift_lhs = F - max_frac;
        constexpr int shift_rhs = Other::frac_bits - max_frac;

        long long lhs_aligned = (shift_lhs == 0) ? raw_ : round_shift(raw_, shift_lhs);
        long long rhs_aligned = (shift_rhs == 0) ? rhs.raw() : round_shift(rhs.raw(), shift_rhs);

        return lhs_aligned == rhs_aligned;
    }

    template<typename Other>
    bool operator!=(const Other& rhs) const {
        return !(*this == rhs);
    }

    // Logarithm operations (input converted to Q16.15, output as Q6.25)
    auto log2() const {
        using Out = FixedPoint<6, 25, Backend>;  // Q6.25 output
        int32_t result = Backend::template log2<total_bits>(raw_, F);
        return Out(result);
    }

    auto logn() const {
        using Out = FixedPoint<6, 25, Backend>;  // Q6.25 output
        int32_t result = Backend::template logn<total_bits>(raw_, F);
        return Out(result);
    }

    auto log10() const {
        using Out = FixedPoint<6, 25, Backend>;  // Q6.25 output
        int32_t result = Backend::template log10<total_bits>(raw_, F);
        return Out(result);
    }

    // Antilogarithm operations (2^x, e^x, 10^x)
    // Input is interpreted as Q6.25, output is Q16.15
    auto antilog2() const {
        using Out = FixedPoint<16, 15, Backend>;  // Q16.15 output
        int32_t result = Backend::template antilog2<total_bits>(raw_, F);
        return Out(result);
    }

    auto antilogn() const {
        using Out = FixedPoint<16, 15, Backend>;  // Q16.15 output
        int32_t result = Backend::template antilogn<total_bits>(raw_, F);
        return Out(result);
    }

    auto antilog10() const {
        using Out = FixedPoint<16, 15, Backend>;  // Q16.15 output
        int32_t result = Backend::template antilog10<total_bits>(raw_, F);
        return Out(result);
    }

    // Power operation: this^exponent
    // Returns a FixedPoint with the same format as this (the base)
    template<typename Other>
    auto pow(const Other& exponent) const {
        using Out = FixedPoint<I, F, Backend>;
        constexpr int Xb = total_bits;
        constexpr int Yb = Other::total_bits;

        using Ax = typename StorageForBits<Xb>::type;
        using By = typename StorageForBits<Yb>::type;

        Ax ax = static_cast<Ax>(raw_);
        By by = static_cast<By>(exponent.raw());
        auto result = Backend::template pow<Xb, Yb>(ax, by, F, Other::frac_bits);
        return Out(result);
    }

    // Square root operation (returns same Q format as input)
    auto sqrt() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template sqrt<total_bits>(raw_, F);
        return Out(result);
    }

    // Reciprocal square root operation: 1/sqrt(x)
    // Returns same Q format as input
    auto rsqrt() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template rsqrt<total_bits>(raw_, F);
        return Out(result);
    }

    // Static array operations (Option 3 API)
    static FixedPoint<I, F, Backend> array_min(const Storage_t<total_bits>* arr, size_t length) {
        auto result = Backend::template array_min<total_bits>(arr, length);
        return FixedPoint<I, F, Backend>(result);
    }

    static FixedPoint<I, F, Backend> array_max(const Storage_t<total_bits>* arr, size_t length) {
        auto result = Backend::template array_max<total_bits>(arr, length);
        return FixedPoint<I, F, Backend>(result);
    }

    // Trigonometric operations (input/output in radians, same Q format)
    auto sin() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template sin<total_bits>(raw_, F);
        return Out(result);
    }

    auto cos() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template cos<total_bits>(raw_, F);
        return Out(result);
    }

    auto tan() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template tan<total_bits>(raw_, F);
        return Out(result);
    }

    auto atan() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template atan<total_bits>(raw_, F);
        return Out(result);
    }

    // Hyperbolic functions (same Q format)
    auto tanh() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template tanh<total_bits>(raw_, F);
        return Out(result);
    }

    // Activation functions (same Q format)
    auto sigmoid() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template sigmoid<total_bits>(raw_, F);
        return Out(result);
    }

    auto relu() const {
        using Out = FixedPoint<I, F, Backend>;
        auto result = Backend::template relu<total_bits>(raw_, F);
        return Out(result);
    }
};

// Short alias
template<int I, int F, typename Backend = ReferenceBackend>
using q = FixedPoint<I,F,Backend>;

// ============================================================================
// FixedPointArray: Wrapper for arrays of fixed-point values (Option 1 API)
// ============================================================================
//
// Provides array-level operations on fixed-point data stored as integer arrays.
// The underlying storage is a pointer to integers (int8_t*, int16_t*, or int32_t*)
// depending on the Q format's total bit width.

template<int I, int F, typename Backend = ReferenceBackend>
class FixedPointArray {
public:
    static constexpr int int_bits = I;
    static constexpr int frac_bits = F;
    static constexpr int total_bits = I + F;
    using Storage = Storage_t<total_bits>;

private:
    Storage* data_;
    size_t length_;

public:
    // Constructor
    FixedPointArray(Storage* data, size_t length)
        : data_(data), length_(length) {}

    // Accessors
    Storage* data() { return data_; }
    const Storage* data() const { return data_; }
    size_t length() const { return length_; }

    // Array element access (returns FixedPoint wrapper)
    FixedPoint<I, F, Backend> operator[](size_t idx) const {
        return FixedPoint<I, F, Backend>(data_[idx]);
    }

    // Array operations
    FixedPoint<I, F, Backend> min() const {
        auto result = Backend::template array_min<total_bits>(data_, length_);
        return FixedPoint<I, F, Backend>(result);
    }

    FixedPoint<I, F, Backend> max() const {
        auto result = Backend::template array_max<total_bits>(data_, length_);
        return FixedPoint<I, F, Backend>(result);
    }

    // In-place array operations
    void shift(int shift_amount) {
        Backend::template array_shift<total_bits>(data_, length_, shift_amount);
    }

    void scale(FixedPoint<I, F, Backend> scale_factor) {
        Backend::template array_scale<total_bits>(data_, length_, scale_factor.raw(), F);
    }

    // Softmax operation (out-of-place, writes to output array)
    void softmax(FixedPointArray<I, F, Backend>& output) const {
        Backend::template softmax<total_bits>(data_, output.data(), length_, F);
    }

    // Vector operations (return scalar FixedPoint results)
    FixedPoint<I, F, Backend> dot_product(const FixedPointArray<I, F, Backend>& other) const {
        auto result = Backend::template dot_product<total_bits>(data_, other.data(), length_, F);
        return FixedPoint<I, F, Backend>(result);
    }

    FixedPoint<I, F, Backend> sum() const {
        auto result = Backend::template array_sum<total_bits>(data_, length_);
        return FixedPoint<I, F, Backend>(result);
    }

    // Element-wise operations (out-of-place, write to output array)
    void elemult(const FixedPointArray<I, F, Backend>& other,
                 FixedPointArray<I, F, Backend>& output) const {
        Backend::template array_elemult<total_bits>(data_, other.data(), output.data(), length_, F);
    }

    void add(const FixedPointArray<I, F, Backend>& other,
             FixedPointArray<I, F, Backend>& output) const {
        Backend::template array_add<total_bits>(data_, other.data(), output.data(), length_);
    }

    void sub(const FixedPointArray<I, F, Backend>& other,
             FixedPointArray<I, F, Backend>& output) const {
        Backend::template array_sub<total_bits>(data_, other.data(), output.data(), length_);
    }

    // Statistical operations (return scalar results)
    FixedPoint<I, F, Backend> mean() const {
        auto result = Backend::template array_mean<total_bits>(data_, length_, F);
        return FixedPoint<I, F, Backend>(result);
    }

    FixedPoint<I, F, Backend> rms() const {
        auto result = Backend::template array_rms<total_bits>(data_, length_, F);
        return FixedPoint<I, F, Backend>(result);
    }

    FixedPoint<I, F, Backend> variance() const {
        auto result = Backend::template array_variance<total_bits>(data_, length_, F);
        return FixedPoint<I, F, Backend>(result);
    }

    FixedPoint<I, F, Backend> stddev() const {
        auto result = Backend::template array_stddev<total_bits>(data_, length_, F);
        return FixedPoint<I, F, Backend>(result);
    }
};

// Short alias for FixedPointArray
template<int I, int F, typename Backend = ReferenceBackend>
using q_array = FixedPointArray<I, F, Backend>;

// ---------- Free helpers (no ambiguous operator overloads) ----------

// Explicit result format helper: fp::mul_as<OUT_I,OUT_F>(a,b)
template<int OUT_I, int OUT_F, typename QA, typename QB>
auto mul_as(const QA& a, const QB& b) {
    return a.template mul<OUT_I, OUT_F>(b);
}

// Explicit result format helper: fp::div_as<OUT_I,OUT_F>(a,b)
template<int OUT_I, int OUT_F, typename QA, typename QB>
auto div_as(const QA& a, const QB& b) {
    return a.template div<OUT_I, OUT_F>(b);
}

// If you *later* want a promoting operator*, define it here,
// but ONLY if you remove the member operator* above to avoid ambiguity.
// (Not included now on purpose.)

} // namespace fp
