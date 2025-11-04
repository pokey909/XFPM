#include <cstdio>
#include <cstdint>

#ifdef __XTENSA__
#include <NatureDSP_Signal.h>
#endif

// Q1.15 helper functions
int16_t q15_from_float(float x) {
    return static_cast<int16_t>(x * 32768.0f);
}

float q15_to_float(int16_t x) {
    return static_cast<float>(x) / 32768.0f;
}

// Reference Q1.15 multiply with proper bit extraction
int16_t q15_mult_reference(int16_t a, int16_t b) {
    int32_t prod = static_cast<int32_t>(a) * static_cast<int32_t>(b);
    // Extract bits [30:15] and round
    int32_t result = (prod + (1 << 14)) >> 15;
    // Saturate to 16-bit range
    if (result > 32767) result = 32767;
    if (result < -32768) result = -32768;
    return static_cast<int16_t>(result);
}

int main() {
    // Test data: simple multiplications
    int16_t arr1[4] = {
        q15_from_float(0.5f),   // 16384
        q15_from_float(0.75f),  // 24576
        q15_from_float(0.25f),  // 8192
        q15_from_float(0.5f)    // 16384
    };

    int16_t arr2[4] = {
        q15_from_float(0.5f),   // 16384
        q15_from_float(0.5f),   // 16384
        q15_from_float(0.5f),   // 16384
        q15_from_float(0.25f)   // 8192
    };

    int16_t output_naturedsp[4];
    int16_t output_reference[4];

    // Expected results:
    // 0.5 * 0.5 = 0.25
    // 0.75 * 0.5 = 0.375
    // 0.25 * 0.5 = 0.125
    // 0.5 * 0.25 = 0.125

    std::printf("Input values:\n");
    for (int i = 0; i < 4; i++) {
        std::printf("  arr1[%d] = %6d (%.6f) * arr2[%d] = %6d (%.6f)\n",
                   i, arr1[i], q15_to_float(arr1[i]),
                   i, arr2[i], q15_to_float(arr2[i]));
    }

    // Compute reference
    for (int i = 0; i < 4; i++) {
        output_reference[i] = q15_mult_reference(arr1[i], arr2[i]);
    }

    std::printf("\nReference implementation results:\n");
    for (int i = 0; i < 4; i++) {
        std::printf("  output[%d] = %6d (%.6f)\n",
                   i, output_reference[i], q15_to_float(output_reference[i]));
    }

#ifdef __XTENSA__
    // Compute using NatureDSP vec_elemult16x16
    vec_elemult16x16(output_naturedsp,
                     const_cast<int16_t*>(arr1),
                     const_cast<int16_t*>(arr2),
                     4);

    std::printf("\nNatureDSP vec_elemult16x16 results:\n");
    for (int i = 0; i < 4; i++) {
        std::printf("  output[%d] = %6d (%.6f)\n",
                   i, output_naturedsp[i], q15_to_float(output_naturedsp[i]));
    }

    std::printf("\nComparison:\n");
    for (int i = 0; i < 4; i++) {
        std::printf("  [%d] Reference: %.6f, NatureDSP: %.6f, Diff: %.6f\n",
                   i,
                   q15_to_float(output_reference[i]),
                   q15_to_float(output_naturedsp[i]),
                   q15_to_float(output_reference[i]) - q15_to_float(output_naturedsp[i]));
    }

    // Also show the raw bit patterns for 0.5 * 0.5
    std::printf("\nDetailed analysis for 0.5 * 0.5:\n");
    int32_t prod = static_cast<int32_t>(arr1[0]) * static_cast<int32_t>(arr2[0]);
    std::printf("  Input: 0x%04X * 0x%04X = 0x%08X (%d)\n",
               static_cast<uint16_t>(arr1[0]),
               static_cast<uint16_t>(arr2[0]),
               static_cast<uint32_t>(prod), prod);
    std::printf("  Bits [31:16] = 0x%04X = %d (%.6f in Q1.15)\n",
               static_cast<uint16_t>(prod >> 16),
               static_cast<int16_t>(prod >> 16),
               q15_to_float(static_cast<int16_t>(prod >> 16)));
    std::printf("  Bits [30:15] = 0x%04X = %d (%.6f in Q1.15)\n",
               static_cast<uint16_t>((prod >> 15) & 0xFFFF),
               static_cast<int16_t>((prod >> 15) & 0xFFFF),
               q15_to_float(static_cast<int16_t>((prod >> 15) & 0xFFFF)));
    std::printf("  NatureDSP gave: 0x%04X = %d (%.6f in Q1.15)\n",
               static_cast<uint16_t>(output_naturedsp[0]),
               output_naturedsp[0],
               q15_to_float(output_naturedsp[0]));
#else
    std::printf("\n(Not compiled for Xtensa - skipping NatureDSP test)\n");
#endif

    return 0;
}
