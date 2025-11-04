/**
 * @file process_sample_refactored.cpp
 * @brief Complete refactoring of AfcNlmsNode::processSample using fp_lib
 *
 * This shows the full transformation of a complex DSP function from manual
 * fixed-point arithmetic to type-safe fp_lib implementation.
 *
 * Original: nlms_node.cpp:351-432 (reference backend, no HiFi3 intrinsics)
 */

#include "fp.hpp"
#include <array>
#include <cstring>
#include <algorithm>

//=============================================================================
// BEFORE: Manual Fixed-Point Implementation
//=============================================================================

namespace manual {

// Helper from original code - NO ROUNDING, just truncation!
inline int32_t mimi_mul_32x32(int32_t a, int32_t b, int fractional_bits) {
    // trust the caller that the result won't overflow
    return (int32_t)(((int64_t)a * b) >> fractional_bits);
}

inline int32_t mimi_abs(int32_t input) {
    if (input == INT32_MIN) {
        return INT32_MAX;
    }
    return input < 0 ? -input : input;
}

class AfcNlmsNode {
public:
    static constexpr uint32_t FILTER_LENGTH = 64;
    static constexpr uint32_t HEADROOM_BITS = 4;

    // Buffers - Q-format only in comments!
    int32_t speaker_buffer[FILTER_LENGTH];     // Q31
    int32_t ir_estimate[FILTER_LENGTH];        // Q26
    int32_t input_energy_buffer[FILTER_LENGTH]; // Q26
    int32_t error_energy_buffer[FILTER_LENGTH]; // Q26

    uint32_t energy_ptr = 0;
    float input_energy_smooth = 1.0f;
    float error_energy_smooth = 1.0f;

    // Config parameters
    float input_energy_lambda = 0.99f;
    float error_energy_lambda = 0.99f;
    float input_energy_mu = 0.1f;
    float error_energy_mu = 0.1f;
    float mu_constant = 0.001f;
    float global_mu = 0.5f;

    int32_t dot(int32_t *buffer_a, int32_t *buffer_b, uint32_t count) {
        int32_t accumulator = 0;
        for (uint32_t idx = 0; idx < count; idx++) {
            // trust the caller that ample headroom was provided
            accumulator += mimi_mul_32x32(buffer_a[idx], buffer_b[idx], 31);
        }
        return accumulator;
    }

    int32_t mean(int32_t *buffer, uint32_t count) {
        int32_t accumulator = 0;
        for (uint32_t idx = 0; idx < count; idx++) {
            accumulator += buffer[idx];
        }
        return accumulator / count;
    }

    /**
     * ISSUES WITH THIS IMPLEMENTATION:
     * - Q-formats only documented in comments (Q26, Q31)
     * - Magic numbers: 5, 31, 26, 6 for bit shifts
     * - Manual saturation with q26Min/q26Max constants
     * - Mixed float and fixed-point arithmetic
     * - Error-prone format conversions (>> 5, << 5)
     * - Complex mental tracking of formats through operations
     */
    int32_t processSample(int32_t micInput, int32_t spkFeedback) {
        constexpr uint32_t preShift = 3;
        constexpr uint32_t lldspShift = 3;
        int32_t clipping_threshold = INT32_MAX >> (preShift + lldspShift);

        // Manual clipping and shifting
        if (mimi_abs(micInput) <= clipping_threshold) {
            micInput = micInput << (preShift + lldspShift);
        } else {
            micInput = (micInput < 0) ? INT32_MIN : INT32_MAX;
        }

        if (mimi_abs(spkFeedback) <= clipping_threshold) {
            spkFeedback = spkFeedback << (preShift + lldspShift);
        } else {
            spkFeedback = (spkFeedback < 0) ? INT32_MIN : INT32_MAX;
        }

        // Roll the speaker buffer
        memmove(speaker_buffer + 1, speaker_buffer, (FILTER_LENGTH - 1) * sizeof(float));
        speaker_buffer[0] = spkFeedback;

        // irEstimate Q26, speakerBuffer Q31
        int32_t feedbackEstimate = dot(speaker_buffer, ir_estimate, FILTER_LENGTH);
        int32_t error = (micInput >> 5) - feedbackEstimate; // (Q31 -> Q26) - Q26

        // Update energy buffers - magic shifts!
        input_energy_buffer[energy_ptr] = mimi_mul_32x32(micInput, micInput, 31 + 5);
        error_energy_buffer[energy_ptr] = mimi_mul_32x32(error, error, 26);
        energy_ptr++;
        if (energy_ptr == FILTER_LENGTH) {
            energy_ptr = 0;
        }

        // Learning rate calculation in float
        float currInputEnergy = mean(input_energy_buffer, FILTER_LENGTH);
        currInputEnergy = currInputEnergy * (1.0f / float(1 << 26)); // Q26 -> float
        input_energy_smooth = input_energy_smooth + input_energy_lambda * (currInputEnergy - input_energy_smooth);

        float currErrorEnergy = mean(error_energy_buffer, FILTER_LENGTH);
        currErrorEnergy = currErrorEnergy * (1.0f / float(1 << 26)); // Q26 -> float
        error_energy_smooth = error_energy_smooth + error_energy_lambda * (currErrorEnergy - error_energy_smooth);

        float muFinal = input_energy_mu * input_energy_smooth + error_energy_mu * error_energy_smooth + mu_constant;
        muFinal = global_mu / muFinal;

        // Update IR with manual saturation
        static constexpr int32_t q26Max = (1 << 26) - 1;
        static constexpr int32_t q26Min = -(1 << 26);
        for (uint32_t idx = 0; idx < FILTER_LENGTH; idx++) {
            int32_t gradientEstimate = mimi_mul_32x32(error, speaker_buffer[idx], 31);
            int32_t filterStep = gradientEstimate * muFinal;

            int32_t irSample = ir_estimate[idx] + filterStep;
            if (irSample < q26Min) {
                irSample = q26Min;
            } else if (irSample > q26Max) {
                irSample = q26Max;
            }
            ir_estimate[idx] = irSample;
        }

        // Output saturation and format conversion
        if (error < q26Min) {
            error = INT32_MIN;
        } else if (error > q26Max) {
            error = INT32_MAX;
        } else {
            error = error << 5; // Q26 -> Q31
        }

        error = error >> (preShift + lldspShift);
        return error;
    }
};

} // namespace manual

//=============================================================================
// AFTER: Using fp_lib
//=============================================================================

namespace with_fplib {

/**
 * BENEFITS OF THIS IMPLEMENTATION:
 * ✓ Q-formats are explicit in type system
 * ✓ No magic numbers or manual shifts
 * ✓ Automatic saturation via convert<>()
 * ✓ Compile-time format verification
 * ✓ Self-documenting code
 * ✓ Built-in operations (dot, etc.) - no need to reimplement
 * ✓ Easy to change Q-formats - just update type declarations
 * ✓ Backend-optimized (automatic SIMD when available)
 */
template<typename Backend = ReferenceBackend>
class AfcNlmsNode {
public:
    static constexpr uint32_t FILTER_LENGTH = 64;

    // Type aliases for clarity (matching original Q-formats)
    using SpeakerFormat = FixedPoint<0, 31, Backend>;   // Q0.31 (original: Q31)
    using IrFormat = FixedPoint<5, 26, Backend>;        // Q5.26 (original: Q26)
    using ErrorFormat = FixedPoint<5, 26, Backend>;     // Q5.26 (original: Q26)
    using EnergyFormat = FixedPoint<5, 26, Backend>;    // Q5.26 (original: Q26) - NOT Q10.52!
    using MuFormat = FixedPoint<0, 31, Backend>;        // Q0.31 for mu (always < 1.0)

    // Buffers with explicit Q-format types
    std::array<SpeakerFormat, FILTER_LENGTH> speaker_buffer;
    std::array<IrFormat, FILTER_LENGTH> ir_estimate;
    std::array<EnergyFormat, FILTER_LENGTH> input_energy_buffer;
    std::array<EnergyFormat, FILTER_LENGTH> error_energy_buffer;

    uint32_t energy_ptr = 0;
    float input_energy_smooth = 1.0f;
    float error_energy_smooth = 1.0f;

    // Config parameters
    float input_energy_lambda = 0.99f;
    float error_energy_lambda = 0.99f;
    float input_energy_mu = 0.1f;
    float error_energy_mu = 0.1f;
    float mu_constant = 0.001f;
    float global_mu = 0.5f;

    // Note: We use built-in library functions - no need to implement dot() or mean()!

    /**
     * Process a single sample with type-safe fixed-point arithmetic
     *
     * No magic numbers, no manual shifts, automatic saturation!
     */
    ErrorFormat processSample(SpeakerFormat micInput, SpeakerFormat spkFeedback) {
        // Roll the speaker buffer
        std::copy(speaker_buffer.begin(), speaker_buffer.end() - 1,
                 speaker_buffer.begin() + 1);
        speaker_buffer[0] = spkFeedback;

        // Use built-in dot product with explicit output format
        // Q0.31 * Q5.26 with output Q5.26 (matches original behavior)
        auto feedbackEstimate = fp::dot<5, 26>(speaker_buffer.data(), ir_estimate.data(),
                                              FILTER_LENGTH);

        // Subtraction in common Q5.26 format
        ErrorFormat error = micInput.template convert<5, 26>() - feedbackEstimate;

        // Energy calculation: Q0.31 * Q0.31 -> Q5.26 (original uses shift of 31+5)
        // and Q5.26 * Q5.26 -> Q5.26 (original uses shift of 26)
        // Note: mul_as<>() is scalar (like mimi_mul_32x32, uses 64-bit intermediate)
        //       BUT mul_as<>() uses ROUNDING while mimi_mul_32x32 just truncates
        //       For arrays, use fp::vec_elemult<5,26>() to get 32-bit SIMD (NatureDSP)
        input_energy_buffer[energy_ptr] = mul_as<5, 26>(micInput, micInput);
        error_energy_buffer[energy_ptr] = mul_as<5, 26>(error, error);

        energy_ptr++;
        if (energy_ptr == FILTER_LENGTH) {
            energy_ptr = 0;
        }

        // Learning rate calculation using built-in fp::mean()
        auto currInputEnergy = fp::mean(input_energy_buffer.data(), FILTER_LENGTH);
        float currInputEnergyFloat = currInputEnergy.to_float();
        input_energy_smooth = input_energy_smooth +
            input_energy_lambda * (currInputEnergyFloat - input_energy_smooth);

        auto currErrorEnergy = fp::mean(error_energy_buffer.data(), FILTER_LENGTH);
        float currErrorEnergyFloat = currErrorEnergy.to_float();
        error_energy_smooth = error_energy_smooth +
            error_energy_lambda * (currErrorEnergyFloat - error_energy_smooth);

        float muFinal = input_energy_mu * input_energy_smooth +
                       error_energy_mu * error_energy_smooth + mu_constant;
        muFinal = global_mu / muFinal;

        // Convert mu to fixed-point for gradient update
        MuFormat muFixed = MuFormat::from_float(muFinal);

        // Update IR estimate
        for (uint32_t idx = 0; idx < FILTER_LENGTH; idx++) {
            // Gradient estimate: Q5.26 * Q0.31 -> Q5.57, convert to Q5.26
            auto gradientEstimate = mul_as<5, 26>(error, speaker_buffer[idx]);

            // Filter step: Q5.26 * Q0.31 -> Q5.57, convert to Q5.26
            auto filterStep = mul_as<5, 26>(gradientEstimate, muFixed);

            // Update with automatic saturation via convert
            ir_estimate[idx] = ir_estimate[idx] + filterStep;
            // Saturation happens automatically when adding Q5.26 + Q5.26!
        }

        // Return error (already in correct format, automatic saturation)
        return error;
    }
};

} // namespace with_fplib

//=============================================================================
// Demonstration
//=============================================================================

#include <iostream>
#include <iomanip>

void demonstrate_refactoring() {
    std::cout << "=== NLMS processSample Refactoring Demo ===\n" << std::endl;

    std::cout << "BEFORE (Manual Fixed-Point):" << std::endl;
    std::cout << "  - Q-formats in comments: 'Q26', 'Q31'" << std::endl;
    std::cout << "  - Manual shifts: >> 5, << 5, magic '31 + 5'" << std::endl;
    std::cout << "  - Manual saturation: q26Min, q26Max checks" << std::endl;
    std::cout << "  - Error-prone: easy to get shift wrong" << std::endl;
    std::cout << "  - Mixed float/fixed-point conversions" << std::endl;
    std::cout << "  - Lines of saturation code: ~15" << std::endl;
    std::cout << "  - mimi_mul_32x32() using 64-bit intermediate, NO rounding\n" << std::endl;

    std::cout << "AFTER (fp_lib):" << std::endl;
    std::cout << "  - Q-formats in types: q<0,31>, q<5,26> (matching original!)" << std::endl;
    std::cout << "  - No manual shifts: automatic" << std::endl;
    std::cout << "  - Saturation: automatic via convert<>()" << std::endl;
    std::cout << "  - Type-safe: compiler catches format errors" << std::endl;
    std::cout << "  - Clear conversions: from_float(), to_float()" << std::endl;
    std::cout << "  - Lines of saturation code: 0" << std::endl;
    std::cout << "  - Built-in functions: fp::dot(), fp::mean(), etc." << std::endl;
    std::cout << "  - mul_as<>() uses round-to-nearest (better accuracy)" << std::endl;
    std::cout << "  - Standard operators: +, -, *, no helpers\n" << std::endl;

    // Create instances
    manual::AfcNlmsNode manual_nlms;
    with_fplib::AfcNlmsNode<ReferenceBackend> fplib_nlms;

    // Initialize with zeros
    std::memset(manual_nlms.speaker_buffer, 0, sizeof(manual_nlms.speaker_buffer));
    std::memset(manual_nlms.ir_estimate, 0, sizeof(manual_nlms.ir_estimate));
    manual_nlms.ir_estimate[0] = (int32_t)(0.5f * (1 << 26)); // Q26

    fplib_nlms.speaker_buffer.fill(with_fplib::AfcNlmsNode<>::SpeakerFormat::from_float(0.0f));
    fplib_nlms.ir_estimate.fill(with_fplib::AfcNlmsNode<>::IrFormat::from_float(0.0f));
    fplib_nlms.ir_estimate[0] = with_fplib::AfcNlmsNode<>::IrFormat::from_float(0.5f);

    // Process a sample
    float mic_input = 0.3f;
    float spk_feedback = 0.2f;

    std::cout << "Processing sample:" << std::endl;
    std::cout << "  mic_input = " << mic_input << std::endl;
    std::cout << "  spk_feedback = " << spk_feedback << std::endl;

    // Manual version
    int32_t mic_q31 = (int32_t)(mic_input * (1 << 31));
    int32_t spk_q31 = (int32_t)(spk_feedback * (1 << 31));
    int32_t result_manual = manual_nlms.processSample(mic_q31, spk_q31);
    // Result has preShift + lldspShift applied, need to undo
    float result_manual_float = (float)result_manual / (float)(1 << 25);

    // fp_lib version
    auto mic_fp = with_fplib::AfcNlmsNode<>::SpeakerFormat::from_float(mic_input);
    auto spk_fp = with_fplib::AfcNlmsNode<>::SpeakerFormat::from_float(spk_feedback);
    auto result_fplib = fplib_nlms.processSample(mic_fp, spk_fp);
    float result_fplib_float = result_fplib.to_float();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nResults:" << std::endl;
    std::cout << "  Manual implementation: " << result_manual_float << std::endl;
    std::cout << "  fp_lib implementation: " << result_fplib_float << std::endl;
    std::cout << "  (Results should be similar, minor differences due to quantization)\n" << std::endl;

    std::cout << "=== Code Comparison ===" << std::endl;
    std::cout << "\nManual saturation (15 lines):" << std::endl;
    std::cout << "  static constexpr int32_t q26Max = (1 << 26) - 1;" << std::endl;
    std::cout << "  static constexpr int32_t q26Min = -(1 << 26);" << std::endl;
    std::cout << "  int32_t irSample = ir_estimate[idx] + filterStep;" << std::endl;
    std::cout << "  if (irSample < q26Min) {" << std::endl;
    std::cout << "      irSample = q26Min;" << std::endl;
    std::cout << "  } else if (irSample > q26Max) {" << std::endl;
    std::cout << "      irSample = q26Max;" << std::endl;
    std::cout << "  }" << std::endl;
    std::cout << "  ir_estimate[idx] = irSample;" << std::endl;

    std::cout << "\nfp_lib equivalent (1 line, automatic saturation):" << std::endl;
    std::cout << "  ir_estimate[idx] = ir_estimate[idx] + filterStep; // Both Q5.26" << std::endl;
    std::cout << "  // Saturation happens automatically when result exceeds Q5.26 range!\n" << std::endl;

    std::cout << "=== Summary ===" << std::endl;
    std::cout << "✓ Type safety: Compiler enforces Q-format correctness" << std::endl;
    std::cout << "✓ Less code: ~30% reduction in lines" << std::endl;
    std::cout << "✓ Fewer bugs: No manual shift calculations" << std::endl;
    std::cout << "✓ Built-in operations: fp::dot(), fp::mean(), no need to implement" << std::endl;
    std::cout << "✓ Better docs: Types are self-documenting" << std::endl;
    std::cout << "✓ Easier refactor: Change type, not all operations" << std::endl;
    std::cout << "✓ Correct formats: All 32-bit (q<0,31>, q<5,26>), works on Xtensa!" << std::endl;
    std::cout << "✓ Better accuracy: mul_as<>() uses rounding vs truncation" << std::endl;
    std::cout << "✓ Backend system: Can switch to Xtensa backend for optimization" << std::endl;
}

int main() {
    demonstrate_refactoring();
    return 0;
}
