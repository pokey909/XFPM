#include <cstdio>
#include "fp.hpp"         // includes all backends
#include <cmath>

namespace test {

static int failures = 0;

inline void expect_near(const char* name, float got, float want, float eps) {
    float err = std::fabs(got - want);
    if (err > eps || std::isnan(got)) {
        std::printf("[FAIL] %s  got=%.9f  want=%.9f  |err|=%.9f  eps=%.9f\n",
                    name, got, want, err, eps);
        ++failures;
    } else {
        std::printf("[ OK ] %s  got=%.9f  want=%.9f  |err|=%.9f\n",
                    name, got, want, err);
    }
}

template<typename QA, typename QB, int OUT_I, int OUT_F>
void check_mul(const char* name, float a, float b, float eps_scale = 1.0f) {
    using OutBackend = typename QA::backend_type;
    using Out        = fp::FixedPoint<OUT_I, OUT_F, OutBackend>;

    auto A = QA::from_float(a);
    auto B = QB::from_float(b);
    auto R = A.template mul<OUT_I, OUT_F>(B);

    float got  = R.to_float();
    float want = a * b;

    const float LSB_out = 1.0f / float(1u << OUT_F);
    const float LSB_a   = 1.0f / float(1u << QA::frac_bits);
    const float LSB_b   = 1.0f / float(1u << QB::frac_bits);

    // If you keep truncation, use 1.0 instead of 0.5 below (or multiply by ~2)
    const float eps = eps_scale * (0.5f*LSB_out + 0.5f*std::fabs(b)*LSB_a + 0.5f*std::fabs(a)*LSB_b);

    test::expect_near(name, got, want, eps);
}

template<typename QA, typename QB, int OUT_I, int OUT_F>
void check_div(const char* name, float a, float b, float eps_scale = 1.0f) {
    using OutBackend = typename QA::backend_type;
    using Out        = fp::FixedPoint<OUT_I, OUT_F, OutBackend>;

    auto A = QA::from_float(a);
    auto B = QB::from_float(b);
    auto R = A.template div<OUT_I, OUT_F>(B);

    float got  = R.to_float();
    float want = a / b;

    const float LSB_out = 1.0f / float(1u << OUT_F);
    const float LSB_a   = 1.0f / float(1u << QA::frac_bits);
    const float LSB_b   = 1.0f / float(1u << QB::frac_bits);

    // Division error propagation is more complex than multiply
    // For now, use a more generous epsilon
    const float eps = eps_scale * (LSB_out + std::fabs(want) * (LSB_a + LSB_b));

    test::expect_near(name, got, want, eps);
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

    test::expect_near(name, got, want, eps);
}

// Quick helper for saturation bounds on Q(I,F) with signed storage
template<int I, int F>
float q_max_value() {
    // storage is signed with total_bits = I+F, so max ≈ (2^(I+F-1)-1)/2^F
    const int total = I + F;
    const int32_t max_raw = (1 << (total - 1)) - 1;
    return static_cast<float>(max_raw) / static_cast<float>(1u << F);
}

} // namespace test

int main() {
    using namespace fp;

    // Select backend based on build target
#ifdef __XTENSA__
    using Backend = XtensaBackend;
#else
    using Backend = ReferenceBackend;
#endif

    // Use selected backend for most tests to exercise tag dispatch + fallbacks
    using q8    = q<1,7,  Backend>;   // total 8  -> int8_t bucket
    using q12   = q<1,11, Backend>;   // total 12 -> int16_t bucket
    using q16   = q<1,15, Backend>;   // total 16 -> int16_t bucket
    using q16b  = q<3,13, Backend>;   // another 16-bit bucket format
    using q32   = q<3,29, Backend>;   // total 32 -> int32_t bucket

    // Also ensure the pure Reference backend works the same
    using q16_ref = q<1,15, ReferenceBackend>;
    using q32_ref = q<3,29, ReferenceBackend>;

    using q1_15 = q<1,15>;  // default backend is ReferenceBackend
    
    // --- 8×8 paths ---
    test::check_mul<q8,   q8,   1,7>("8x8->8  (Backend prio-2)",  0.50f, 0.75f);
    test::check_mul<q8,   q8,   1,15>("8x8->16 (fallback)",             0.50f, 0.75f);

    // --- 16×16 paths ---
    test::check_mul<q16,  q16,  1,15>("16x16->16 (Backend prio-1)", 0.50f, 0.75f);
    test::check_mul<q16,  q16,  3,29>("16x16->32 (fallback)",             0.50f, 0.75f);

    // --- 12-bit total (maps to 16-bit bucket internally) ---
    test::check_mul<q12,  q12,  3,13>("12b×12b -> 16b bucket", 0.30f, 0.40f);

    // --- Mixed-width operands ---
    test::check_mul<q8,   q16,  3,13>("8b×16b -> 16b",          0.40f, 0.60f);
    test::check_mul<q16,  q8,   3,13>("16b×8b -> 16b",          0.40f, 0.60f);
    test::check_mul<q16,  q32,  3,29>("16b×32b -> 32b",         0.40f, 0.60f);

    // --- Different OUT_F (forces nontrivial right shift) ---
    // frac_in = 15+15 = 30; OUT_F = 10 → shift = 20
    test::check_mul<q16,  q16,  3,10>("16x16 -> Q3.10 (shift test)", 0.50f, 0.25f, 2.0f);

    // --- Reference backend parity ---
    test::check_mul<q16_ref, q16_ref, 1,15>("Reference 16x16->16", 0.50f, 0.75f);
    test::check_mul<q16_ref, q16_ref, 3,29>("Reference 16x16->32", 0.50f, 0.75f);

    // --- Saturation sanity (values near representable max) ---
    {
        // Q1.15 max is just below 1.0
        float q115_max = test::q_max_value<1,15>();
        auto A = q16::from_float(q115_max);
        auto B = q16::from_float(q115_max);
        auto R = A.template mul<1,15>(B);  // stays below max; numerical <= max
        float got = R.to_float();
        // Expected is q115_max^2 (~<1), within a couple LSBs
        float want = q115_max * q115_max;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("Saturation-safe multiply near max (Q1.15)",
                          got, want, 2.0f * lsb);
    }

    // --- Mixed Q formats with different integer spans (Q3.13 * Q1.15 -> Q3.13) ---
    test::check_mul<q16b, q16,  3,13>("Q3.13 * Q1.15 -> Q3.13", 0.80f, 0.50f);

    // --- A small negative case to confirm signed behavior ---
    test::check_mul<q16,  q16,  1,15>("Signed: (-0.5)*0.75",  -0.50f, 0.75f);
    test::check_mul<q8,   q16,  3,13>("Signed: (-0.3)*0.6",   -0.30f, 0.60f);

    q1_15 a = q1_15::from_float(0.5f);
    q1_15 b = q1_15::from_float(0.75f);

    auto same = a * b;                   // Q1.15 result
    auto hi   = fp::mul_as<3,29>(a,b);   // Q3.29 result

    printf("%f  %f\n", same.to_float(), hi.to_float());

    // --- Division tests ---
    std::puts("\n--- Division Tests ---");

    // Basic division tests (results must fit in output format)
    // Q1.15 can represent [-1.0, 1.0), so 0.375/0.50 = 0.75 fits
    test::check_div<q16,  q16,  1,15>("16÷16->16 (Backend prio-1)", 0.375f, 0.50f);
    // Q1.7 can represent [-1.0, 1.0), so 0.25/0.50 = 0.5 fits
    test::check_div<q8,   q8,   1,7> ("8÷8->8  (Backend prio-2)",  0.25f, 0.50f);
    // Q3.29 can represent up to ~4.0, so 0.75/0.50 = 1.5 fits
    test::check_div<q16,  q16,  3,29>("16÷16->32 (fallback)",            0.75f, 0.50f);

    // Mixed width operands (results fit in Q1.15 and Q1.7)
    test::check_div<q16,  q8,   1,15>("16÷8->16",  0.30f, 0.60f);  // 0.5
    test::check_div<q8,   q16,  1,7> ("8÷16->8",   0.20f, 0.40f);  // 0.5

    // Negative values (results fit in Q1.15 range)
    test::check_div<q16,  q16,  1,15>("Signed: (-0.375)/0.5",  -0.375f, 0.50f);  // -0.75
    test::check_div<q16,  q16,  1,15>("Signed: 0.375/(-0.5)",   0.375f, -0.50f);  // -0.75
    test::check_div<q16,  q16,  1,15>("Signed: (-0.375)/(-0.5)", -0.375f, -0.50f);  // 0.75

    // Test operator/ with values that fit
    {
        auto c = q16::from_float(0.375f);
        auto d = q16::from_float(0.50f);
        auto quotient = c / d;  // Should use operator/
        float got = quotient.to_float();
        float want = 0.375f / 0.50f;  // = 0.75
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("Operator / test", got, want, 2.0f * lsb);
    }

    // Test free function div_as with result that needs more integer bits
    {
        auto e = q16::from_float(0.75f);
        auto f = q16::from_float(0.50f);
        auto quotient32 = fp::div_as<3,29>(e, f);  // Q3.29 can hold 1.5
        float got = quotient32.to_float();
        float want = 0.75f / 0.50f;  // = 1.5
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        test::expect_near("div_as<3,29> test", got, want, 2.0f * lsb);
    }

    // Reference backend parity
    test::check_div<q16_ref, q16_ref, 1,15>("Reference 16÷16->16", 0.375f, 0.50f);
    test::check_div<q32_ref, q32_ref, 3,29>("Reference 32÷32->32", 0.80f, 0.40f);

    // --- Logarithm tests ---
    std::puts("\n--- Logarithm Tests ---");

    // Test log2()
    {
        auto x = q16::from_float(0.5f);  // log2(0.5) = -1.0
        auto result = x.log2();
        float got = result.to_float();
        float want = std::log2(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);  // Q6.25
        test::expect_near("log2(0.5)", got, want, 10.0f * lsb);
    }

    {
        auto x = q16::from_float(0.25f);  // log2(0.25) = -2.0
        auto result = x.log2();
        float got = result.to_float();
        float want = std::log2(0.25f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        test::expect_near("log2(0.25)", got, want, 10.0f * lsb);
    }

    // Test logn() (natural log)
    {
        auto x = q16::from_float(0.5f);  // ln(0.5) ≈ -0.693
        auto result = x.logn();
        float got = result.to_float();
        float want = std::log(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        test::expect_near("logn(0.5)", got, want, 10.0f * lsb);
    }

    // Test log10()
    {
        auto x = q16::from_float(0.1f);  // log10(0.1) = -1.0
        auto result = x.log10();
        float got = result.to_float();
        float want = std::log10(0.1f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        test::expect_near("log10(0.1)", got, want, 1000.0f * lsb);
    }

    // Test with 32-bit input
    {
        auto x = q32::from_float(0.5f);
        auto result = x.log2();
        float got = result.to_float();
        float want = std::log2(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 25);
        test::expect_near("log2(0.5) with Q32", got, want, 10.0f * lsb);
    }

    // Test error handling for zero/negative
    {
        auto x = q16::from_float(-0.5f);  // negative input
        auto result = x.log2();
        // Result should be most negative value (0x80000000)
        // which as Q6.25 is -64.0
        float got = result.to_float();
        float want = -64.0f;  // Approximate expected error value
        test::expect_near("log2(-0.5) error handling", got, want, 1.0f);
    }

    // --- Antilogarithm tests ---
    std::puts("\n--- Antilogarithm Tests ---");

    // Test antilog2 (2^x)
    // Input is interpreted as Q6.25, output is Q16.15
    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(1.0f);  // 2^1 = 2
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 2.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog2(1.0) = 2.0", got, want, 2.0f * lsb);
    }

    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(2.0f);  // 2^2 = 4
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 4.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog2(2.0) = 4.0", got, want, 2.0f * lsb);
    }

    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(-1.0f);  // 2^-1 = 0.5
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 0.5f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog2(-1.0) = 0.5", got, want, 2.0f * lsb);
    }

    // Test antilogn (e^x)
    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(1.0f);  // e^1 ≈ 2.718
        auto result = x.antilogn();
        float got = result.to_float();
        float want = std::exp(1.0f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilogn(1.0) ≈ 2.718", got, want, 2.0f * lsb);
    }

    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(0.0f);  // e^0 = 1
        auto result = x.antilogn();
        float got = result.to_float();
        float want = 1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilogn(0.0) = 1.0", got, want, 2.0f * lsb);
    }

    // Test antilog10 (10^x)
    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(1.0f);  // 10^1 = 10
        auto result = x.antilog10();
        float got = result.to_float();
        float want = 10.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog10(1.0) = 10.0", got, want, 2.0f * lsb);
    }

    {
        using q6_25 = q<6,25, Backend>;
        auto x = q6_25::from_float(0.0f);  // 10^0 = 1
        auto result = x.antilog10();
        float got = result.to_float();
        float want = 1.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog10(0.0) = 1.0", got, want, 2.0f * lsb);
    }

    // Test with 32-bit input (priority 1 path for Xtensa)
    {
        using q6_25_32 = q<6,25, Backend>;
        auto x = q6_25_32::from_float(3.0f);  // 2^3 = 8
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog2(3.0) with Q32 (Xtensa prio-1)", got, want, 2.0f * lsb);
    }

    // Test Reference backend
    {
        using q6_25_ref = q<6,25, ReferenceBackend>;
        auto x = q6_25_ref::from_float(2.0f);  // 2^2 = 4
        auto result = x.antilog2();
        float got = result.to_float();
        float want = 4.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("antilog2(2.0) Reference backend", got, want, 2.0f * lsb);
    }

    // --- Square Root tests ---
    std::puts("\n--- Square Root Tests ---");

    // Basic sqrt tests with common values
    {
        auto x = q16::from_float(0.25f);  // sqrt(0.25) = 0.5
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.25f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("sqrt(0.25)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16::from_float(0.5f);  // sqrt(0.5) ≈ 0.707
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("sqrt(0.5)", got, want, 2.0f * lsb);
    }

    {
        auto x = q16::from_float(0.16f);  // sqrt(0.16) = 0.4
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.16f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("sqrt(0.16)", got, want, 2.0f * lsb);
    }

    // Test with Q3.29 format (32-bit, more precision)
    {
        auto x = q32::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.25f);
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        test::expect_near("sqrt(0.25) with Q32", got, want, 2.0f * lsb);
    }

    {
        auto x = q32::from_float(2.0f);  // sqrt(2) ≈ 1.414
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(2.0f);
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        test::expect_near("sqrt(2.0) with Q32", got, want, 2.0f * lsb);
    }

    // Test with Q1.7 format (8-bit)
    {
        auto x = q8::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.25f);
        const float lsb = 1.0f / static_cast<float>(1u << 7);
        test::expect_near("sqrt(0.25) with Q8", got, want, 2.0f * lsb);
    }

    // Test with zero
    {
        auto x = q16::from_float(0.0f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = 0.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("sqrt(0.0)", got, want, lsb);
    }

    // Test error handling for negative input
    {
        auto x = q16::from_float(-0.25f);
        auto result = x.sqrt();
        // Result should be most negative value for int16_t (0x8000)
        // which as Q1.15 is -1.0
        float got = result.to_float();
        float want = -1.0f;
        test::expect_near("sqrt(-0.25) error handling", got, want, 0.1f);
    }

    // Test error handling for negative input with Q32
    {
        auto x = q32::from_float(-1.0f);
        auto result = x.sqrt();
        // Result should be most negative value for int32_t (0x80000000)
        // which as Q3.29 is approximately -4.0
        float got = result.to_float();
        float want = -4.0f;
        test::expect_near("sqrt(-1.0) error handling Q32", got, want, 0.1f);
    }

    // Reference backend parity
    {
        auto x = q16_ref::from_float(0.25f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.25f);
        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("Reference sqrt(0.25)", got, want, 2.0f * lsb);
    }

    {
        auto x = q32_ref::from_float(0.5f);
        auto result = x.sqrt();
        float got = result.to_float();
        float want = std::sqrt(0.5f);
        const float lsb = 1.0f / static_cast<float>(1u << 29);
        test::expect_near("Reference sqrt(0.5)", got, want, 2.0f * lsb);
    }

    // --- Power tests ---
    std::puts("\n--- Power Tests ---");

    // Test case 1: pow(2.0, 3.0) = 8.0
    // Need Q format that can hold 8.0 for result, and 3.0 for exponent
    // Use Q4.12 for base (can hold up to ~16) and Q3.13 for exponent (can hold up to ~4)
    {
        using q4_12 = q<4,12, Backend>;
        using q3_13 = q<3,13, Backend>;  // Can hold [-4.0, 4.0)
        test::check_pow<q4_12, q3_13>("pow(2.0, 3.0) = 8.0", 2.0f, 3.0f, 2.0f);
    }

    // Test case 2: pow(0.5, 2.0) = 0.25
    // Use Q3.13 for exponent to hold 2.0
    {
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q16, q3_13>("pow(0.5, 2.0) = 0.25", 0.5f, 2.0f, 2.0f);
    }

    // Test case 3: pow(4.0, 0.5) = 2.0 (square root of 4)
    {
        using q4_12 = q<4,12, Backend>;
        test::check_pow<q4_12, q16>("pow(4.0, 0.5) = 2.0", 4.0f, 0.5f, 2.0f);
    }

    // Test case 4: pow(0.25, 2.0) = 0.0625
    {
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q16, q3_13>("pow(0.25, 2.0) = 0.0625", 0.25f, 2.0f, 2.0f);
    }

    // Test case 5: pow(3.0, 2.0) = 9.0
    // Result 9.0 needs Q5.11 (can hold [-16.0, 16.0))
    {
        using q5_11 = q<5,11, Backend>;
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q5_11, q3_13>("pow(3.0, 2.0) = 9.0", 3.0f, 2.0f, 2.0f);
    }

    // Test case 6: pow(0.5, 0.5) ≈ 0.707
    test::check_pow<q16, q16>("pow(0.5, 0.5) ≈ 0.707", 0.5f, 0.5f, 2.0f);

    // Test case 7: Fractional exponent - pow(8.0, 1/3) ≈ 2.0 (cube root)
    {
        using q4_12 = q<4,12, Backend>;
        float one_third = 1.0f / 3.0f;
        test::check_pow<q4_12, q16>("pow(8.0, 1/3) ≈ 2.0", 8.0f, one_third, 2.0f);
    }

    // Test case 8: pow(2.0, 0.0) = 1.0 (anything to power 0)
    // Base 2.0 doesn't fit in Q1.15, use Q3.13 for base
    {
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q3_13, q16>("pow(2.0, 0.0) = 1.0", 2.0f, 0.0f, 2.0f);
    }

    // Test case 9: pow(1.0, x) = 1.0 for any x (but x must fit)
    // Exponent 3.0 needs Q3.13
    {
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q16, q3_13>("pow(1.0, 3.0) = 1.0", 1.0f, 3.0f, 2.0f);
    }

    // Test case 10: Negative exponent - pow(2.0, -1.0) = 0.5
    // Base 2.0 needs Q3.13, exponent -1.0 fits in Q1.15
    {
        using q3_13 = q<3,13, Backend>;
        test::check_pow<q3_13, q16>("pow(2.0, -1.0) = 0.5", 2.0f, -1.0f, 2.0f);
    }

    // Test zero/negative base handling
    {
        using q3_13 = q<3,13, Backend>;
        auto base_zero = q16::from_float(0.0f);
        auto exp = q3_13::from_float(2.0f);  // Use Q3.13 for exponent 2.0
        auto result = base_zero.pow(exp);
        float got = result.to_float();
        test::expect_near("pow(0.0, 2.0) = 0.0 (error handling)", got, 0.0f, 0.001f);
    }

    {
        using q3_13 = q<3,13, Backend>;
        auto base_neg = q16::from_float(-0.5f);
        auto exp = q3_13::from_float(2.0f);  // Use Q3.13 for exponent 2.0
        auto result = base_neg.pow(exp);
        float got = result.to_float();
        test::expect_near("pow(-0.5, 2.0) = 0.0 (error handling)", got, 0.0f, 0.001f);
    }

    // Test with 32-bit base (priority 1 path)
    // Result 8.0 needs Q5.27 (can hold [-16.0, 16.0))
    {
        using q5_27 = q<5,27, Backend>;
        using q3_13 = q<3,13, Backend>;
        auto base = q5_27::from_float(2.0f);
        auto exp = q3_13::from_float(3.0f);  // Use Q3.13 for exponent 3.0
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 27);
        test::expect_near("pow(2.0, 3.0) with Q32 base (Xtensa prio-1)", got, want, 2.0f * lsb);
    }

    // Test with 8-bit base (fallback to priority 0)
    {
        using q3_13 = q<3,13, Backend>;
        auto base = q8::from_float(0.5f);
        auto exp = q3_13::from_float(2.0f);
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 0.25f;
        const float lsb = 1.0f / static_cast<float>(1u << 7);
        test::expect_near("pow(0.5, 2.0) with Q8 (fallback)", got, want, 2.0f * lsb);
    }

    // Test Reference backend parity
    // Base 2.0 and result 8.0 need Q4.12 (can hold [-8.0, 8.0))
    {
        using q4_12_ref = q<4,12, ReferenceBackend>;
        using q3_13_ref = q<3,13, ReferenceBackend>;
        auto base = q4_12_ref::from_float(2.0f);
        auto exp = q3_13_ref::from_float(3.0f);
        auto result = base.pow(exp);
        float got = result.to_float();
        float want = 8.0f;
        const float lsb = 1.0f / static_cast<float>(1u << 12);
        test::expect_near("pow(2.0, 3.0) Reference backend", got, want, 10.0f * lsb);
    }

    // --- Array Min/Max tests ---
    std::puts("\n--- Array Min/Max Tests ---");

    // Test Option 1 API: FixedPointArray
    {
        int16_t data[] = {
            q16::from_float(0.5f).raw(),
            q16::from_float(-0.25f).raw(),
            q16::from_float(0.75f).raw(),
            q16::from_float(0.1f).raw(),
            q16::from_float(-0.5f).raw()
        };

        fp::FixedPointArray<1, 15, Backend> arr(data, 5);

        auto min_val = arr.min();
        auto max_val = arr.max();

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.5f;
        float want_max = 0.75f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("FixedPointArray::min()", got_min, want_min, 2.0f * lsb);
        test::expect_near("FixedPointArray::max()", got_max, want_max, 2.0f * lsb);
    }

    // Test Option 3 API: Static methods on FixedPoint
    {
        int16_t data[] = {
            q16::from_float(0.3f).raw(),
            q16::from_float(0.9f).raw(),
            q16::from_float(-0.2f).raw(),
            q16::from_float(0.6f).raw()
        };

        auto min_val = q16::array_min(data, 4);
        auto max_val = q16::array_max(data, 4);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.2f;
        float want_max = 0.9f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("q16::array_min() static", got_min, want_min, 2.0f * lsb);
        test::expect_near("q16::array_max() static", got_max, want_max, 2.0f * lsb);
    }

    // Test with 8-bit arrays (priority 2 path for Xtensa)
    {
        using q8 = q<1, 7, Backend>;
        int8_t data[] = {
            q8::from_float(0.5f).raw(),
            q8::from_float(0.25f).raw(),
            q8::from_float(-0.5f).raw(),
            q8::from_float(0.75f).raw()
        };

        auto min_val = q8::array_min(data, 4);
        auto max_val = q8::array_max(data, 4);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.5f;
        float want_max = 0.75f;

        const float lsb = 1.0f / static_cast<float>(1u << 7);
        test::expect_near("q8 array_min (Xtensa prio-2)", got_min, want_min, 2.0f * lsb);
        test::expect_near("q8 array_max (Xtensa prio-2)", got_max, want_max, 2.0f * lsb);
    }

    // Test with Reference backend
    {
        using q16_ref = q<1, 15, ReferenceBackend>;
        int16_t data[] = {
            q16_ref::from_float(0.2f).raw(),
            q16_ref::from_float(-0.7f).raw(),
            q16_ref::from_float(0.4f).raw()
        };

        auto min_val = q16_ref::array_min(data, 3);
        auto max_val = q16_ref::array_max(data, 3);

        float got_min = min_val.to_float();
        float got_max = max_val.to_float();
        float want_min = -0.7f;
        float want_max = 0.4f;

        const float lsb = 1.0f / static_cast<float>(1u << 15);
        test::expect_near("Reference backend array_min", got_min, want_min, 2.0f * lsb);
        test::expect_near("Reference backend array_max", got_max, want_max, 2.0f * lsb);
    }

    // Summary
    if (test::failures == 0) {
        std::puts("\n==== ALL TESTS PASSED ====");
        return 0;
    } else {
        std::printf("\n==== %d TEST(S) FAILED ====\n", test::failures);
        return 1;
    }
}
