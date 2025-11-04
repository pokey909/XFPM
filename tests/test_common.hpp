#pragma once
#include "../fp.hpp"
#include <cstdio>
#include <cmath>

namespace fp {
namespace test {

// Conditional backend selection based on build target
#ifdef __XTENSA__
using Backend = XtensaBackend;
#else
using Backend = ReferenceBackend;
#endif

// Test state
inline int failures = 0;

// Helper: Expect two values to be near each other
inline void expect_near(const char* name, float got, float want, float eps) {
    float err = std::fabs(got - want);
    if (err <= eps) {
        std::printf("[ OK ] %s  got=%.9f  want=%.9f  |err|=%.9f\n", name, got, want, err);
    } else {
        std::printf("[FAIL] %s  got=%.9f  want=%.9f  |err|=%.9f  eps=%.9f\n",
                    name, got, want, err, eps);
        failures++;
    }
}

// Template helpers for multiply and divide tests
template<typename QA, typename QB, int OutI, int OutF>
void check_mul(const char* name, float a, float b, float eps_scale = 2.0f) {
    auto A = QA::from_float(a);
    auto B = QB::from_float(b);
    auto R = fp::mul_as<OutI, OutF>(A, B);

    float got  = R.to_float();
    float want = a * b;

    const float LSB_a   = 1.0f / float(1u << QA::frac_bits);
    const float LSB_b   = 1.0f / float(1u << QB::frac_bits);
    const float LSB_out = 1.0f / float(1u << OutF);

    // Basic error bound from quantization with scaling factor
    const float eps = eps_scale * (LSB_out + std::fabs(want) * (LSB_a + LSB_b));

    expect_near(name, got, want, eps);
}

template<typename QA, typename QB, int OutI, int OutF>
void check_div(const char* name, float a, float b, float eps_scale = 2.0f) {
    auto A = QA::from_float(a);
    auto B = QB::from_float(b);
    auto R = fp::div_as<OutI, OutF>(A, B);

    float got  = R.to_float();
    float want = a / b;

    const float LSB_a   = 1.0f / float(1u << QA::frac_bits);
    const float LSB_b   = 1.0f / float(1u << QB::frac_bits);
    const float LSB_out = 1.0f / float(1u << OutF);

    // Division error is more complex, use generous bound
    const float eps = eps_scale * (LSB_out + std::fabs(want) * (LSB_a + LSB_b));

    expect_near(name, got, want, eps);
}

template<typename QA, typename QB>
void check_pow(const char* name, float base, float exponent, float eps_scale = 1.0f) {
    auto A = QA::from_float(base);
    auto B = QB::from_float(exponent);
    auto R = A.pow(B);

    float got  = R.to_float();
    float want = std::pow(base, exponent);

    const float LSB_base = 1.0f / float(1u << QA::frac_bits);
    const float LSB_exp  = 1.0f / float(1u << QB::frac_bits);

    // Power function error propagation is complex, use generous epsilon
    const float eps = eps_scale * (LSB_base + std::fabs(want) * (LSB_base + LSB_exp));

    expect_near(name, got, want, eps);
}

// Quick helper for saturation bounds on Q(I,F) with signed storage
template<int I, int F>
float sat_max() {
    return static_cast<float>((1 << I) - 1) + static_cast<float>((1u << F) - 1) / static_cast<float>(1u << F);
}

} // namespace test
} // namespace fp
