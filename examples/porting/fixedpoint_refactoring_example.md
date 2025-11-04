# Fixed-Point Library Refactoring Example

This document demonstrates how to refactor manual fixed-point arithmetic code to use the `fp_lib` library, using the NLMS (Normalized Least Mean Squares) algorithm as a real-world example.

## Overview

The NLMS implementation in `nlms_node.cpp` contains manual fixed-point arithmetic with:
- Explicit Q-format tracking in comments (Q26, Q31, etc.)
- Manual shift operations for format conversions
- Manual saturation logic
- Helper functions like `mimi_mul_32x32` for fixed-point multiplication
- Manual implementations of common operations (dot product, mean, etc.)
- Complex bit-width management

This refactoring shows how the `fp_lib` library eliminates these pain points through five detailed examples.

---

## Example 1: Dot Product Operation

### Before: Manual Fixed-Point (Reference Backend)

```cpp
// From nlms_node.cpp:320-327
int32_t AfcNlmsNode::dot(int32_t *buffer_a, int32_t *buffer_b, uint32_t count) {
    int32_t accumulator = 0;
    for (uint32_t idx = 0; idx < count; idx++) {
        // trust the caller that ample headroom was provided
        accumulator += mimi_mul_32x32(buffer_a[idx], buffer_b[idx], 31);
    }
    return accumulator;
}
```

**Issues:**
- Q-format is implicit (caller must know buffer_a is Q31, buffer_b is Q26)
- `mimi_mul_32x32` with magic number `31` for fractional bits
- No compile-time type safety
- Saturation/overflow handling is manual and error-prone
- Comment indicates trust in caller for headroom

### After: Using fp_lib

```cpp
#include "fp.hpp"

// The library already provides dot product!
// Just use it directly with your fixed-point arrays

// Usage example (matching original NLMS Q-formats):
q<0,31> speaker_buffer[256];  // Q0.31 format (original: Q31)
q<5,26> ir_estimate[256];     // Q5.26 format (original: Q26)

// Specify output format explicitly to match original behavior:
// Original does Q31 * Q26 with shift of 31 -> Q26 output
auto result = fp::dot<5, 26>(speaker_buffer, ir_estimate, 256);
// Result is q<5,26> - stays in 32-bit format, works on Xtensa!

// Without explicit format, result would be q<5,57> (needs 62 bits - NOT supported!)
// auto result_auto = fp::dot(speaker_buffer, ir_estimate, 256); // Don't do this!
```

**Benefits:**
- ✅ Built-in function - no need to implement yourself!
- ✅ Q-format is explicit in type signatures
- ✅ Compile-time verification of format compatibility
- ✅ Automatic saturation and rounding
- ✅ No magic numbers or manual shift calculations
- ✅ Self-documenting code
- ✅ Backend-optimized (uses SIMD when available)

---

## Example 2: Energy Calculation with Format Conversion

### Before: Manual Fixed-Point

```cpp
// From nlms_node.cpp:385-386
// Update energy buffers used for normalization
// inputEnergyBuffer Q26
inputEnergyBuffer.i32[energyPtr] = mimi_mul_32x32(micInput, micInput, 31 + 5);
errorEnergyBuffer.i32[energyPtr] = mimi_mul_32x32(error, error, 26);
```

**Issues:**
- Magic shift values: `31 + 5` and `26`
- Q-format only documented in comments
- Easy to get shift calculation wrong
- Format conversion `Q31 * Q31 -> Q26` requires mental math

### After: Using fp_lib

```cpp
#include "fp.hpp"

// Explicit format conversion using mul_as<>
q<0,31> micInput;     // Q0.31 format (original: Q31)
q<5,26> error;        // Q5.26 format (original: Q26)

// Energy calculation with explicit output format matching original Q26
// Q0.31 * Q0.31 -> Q5.26 (requires shift of 31+5, handled by mul_as)
inputEnergyBuffer[energyPtr] = mul_as<5,26>(micInput, micInput);  // Q5.26
errorEnergyBuffer[energyPtr] = mul_as<5,26>(error, error);        // Q5.26

// Note: Output stays in 32-bit Q5.26 format, NOT Q10.52 (which would need 62 bits!)
// This matches the original code and works on Xtensa hardware
```

**Benefits:**
- ✅ Q-format is part of the type system
- ✅ No manual shift calculation
- ✅ Explicit output format with `mul_as<I,F>()`
- ✅ Automatic precision handling (64-bit intermediate like `mimi_mul_32x32`)

**Note**: This example shows **scalar** operations (sample-by-sample). For **array** operations, use `fp::vec_elemult<5,26>(result, array_a, array_b, count)` to leverage 32-bit SIMD on Xtensa.

---

## Example 3: Complex Processing with Multiple Format Conversions

### Before: Manual Fixed-Point (from processSample)

```cpp
// From nlms_node.cpp:379-390
// irEstimate Q26, speakerBuffer Q31
int32_t feedbackEstimate = dot(speakerBuffer.i32, irEstimate.i32, config.filterLength);
int32_t error = (micInput >> 5) - feedbackEstimate; // (Q31 -> Q26) - Q26

// update energy buffers used for normalization
// inputEnergyBuffer Q26
inputEnergyBuffer.i32[energyPtr] = mimi_mul_32x32(micInput, micInput, 31 + 5);
errorEnergyBuffer.i32[energyPtr] = mimi_mul_32x32(error, error, 26); // Q26 * Q26 -> Q26
energyPtr++;
if(energyPtr == config.filterLength) {
    energyPtr = 0;
}
```

**Issues:**
- Manual right shift `>> 5` for format conversion
- Complex comment tracking: `(Q31 -> Q26) - Q26`
- Shift value `31 + 5` requires understanding the entire format chain
- Buffer types don't enforce Q-format

### After: Using fp_lib

```cpp
#include "fp.hpp"

// Type-safe buffers with explicit Q-format (matching original code)
std::array<q<0,31>, 256> speakerBuffer;   // Q0.31 (original: Q31)
std::array<q<5,26>, 256> irEstimate;      // Q5.26 (original: Q26)
std::array<q<5,26>, 256> inputEnergyBuffer;  // Q5.26 (original: Q26)
std::array<q<5,26>, 256> errorEnergyBuffer;  // Q5.26 (original: Q26)

// Processing with automatic format handling
q<0,31> micInput;  // Q0.31 (original: Q31)

// Use built-in dot product with explicit output format
// Q0.31 * Q5.26 -> Q5.26 (matches original shift behavior)
auto feedbackEstimate = fp::dot<5, 26>(speakerBuffer.data(), irEstimate.data(),
                                       config.filterLength);

// Error calculation - library handles format matching
q<5,26> error = micInput.template convert<5,26>() - feedbackEstimate;

// Energy calculation with explicit Q5.26 output (fits in 32 bits)
inputEnergyBuffer[energyPtr] = mul_as<5,26>(micInput, micInput);
errorEnergyBuffer[energyPtr] = mul_as<5,26>(error, error);

energyPtr++;
if(energyPtr == config.filterLength) {
    energyPtr = 0;
}
```

**Benefits:**
- ✅ No manual bit shifts
- ✅ Format conversions are explicit and type-safe
- ✅ Compile-time errors for incompatible operations
- ✅ Self-documenting types

---

## Example 4: Mean Calculation (Built-in Function)

### Before: Manual Fixed-Point

```cpp
// From nlms_node.cpp:329-336
int32_t AfcNlmsNode::mean(int32_t *buffer, uint32_t count) {
    int32_t accumulator = 0.0f;
    for (uint32_t idx = 0; idx < count; idx++) {
        accumulator += buffer[idx];
    }
    return accumulator / count;
}

// Usage in learning rate calculation (lines 393-400)
float currInputEnergy = mean(inputEnergyBuffer.i32, config.filterLength);
currInputEnergy = currInputEnergy * (1.0f / float(1 << 26)); // Q26 -> float
inputEnergySmooth = inputEnergySmooth + config.inputEnergyLambda * (currInputEnergy - inputEnergySmooth);
```

**Issues:**
- Need to manually implement mean function
- Manual Q-format conversion: `* (1.0f / float(1 << 26))`
- Separate function needed for each Q-format
- Magic number `26` tied to specific format

### After: Using fp_lib

```cpp
// No need to implement mean - it's built-in!

// Usage with built-in fp::mean()
auto currInputEnergy = fp::mean(inputEnergyBuffer.data(), config.filterLength);
float currInputEnergyFloat = currInputEnergy.to_float();  // Type-safe conversion
inputEnergySmooth = inputEnergySmooth + config.inputEnergyLambda * (currInputEnergyFloat - inputEnergySmooth);
```

**Benefits:**
- ✅ Built-in function - no manual implementation needed
- ✅ Works with any Q-format automatically
- ✅ Type-safe `.to_float()` conversion
- ✅ No magic numbers
- ✅ Backend-optimized (SIMD when available)

---

## Example 5: Gradient Update with Saturation

### Before: Manual Fixed-Point

```cpp
// From nlms_node.cpp:405-419
// update the IR according to the gradient estimate
static constexpr int32_t q26Max = (1 << 26) - 1;
static constexpr int32_t q26Min = -(1 << 26);
for (uint32_t idx = 0; idx < config.filterLength; idx++) {
    int32_t gradientEstimate = mimi_mul_32x32(error, speakerBuffer.i32[idx], 31); // Q26 * Q31 -> Q26
    int32_t filterStep = gradientEstimate * muFinal;  // muFinal is float!

    int32_t irSample = irEstimate.i32[idx] + filterStep;
    if (irSample < q26Min) {
        irSample = q26Min;
    } else if (irSample > q26Max) {
        irSample = q26Max;
    }
    irEstimate.i32[idx] = irSample;
}
```

**Issues:**
- Manual saturation constants `q26Max`, `q26Min`
- Mixed float and fixed-point (`muFinal` is float)
- Easy to forget saturation checks
- Magic number `31` in multiply
- Manual min/max clamping

### After: Using fp_lib

```cpp
#include "fp.hpp"

// Strongly-typed fixed-point values (matching original Q-formats)
q<5,26> error;                          // Q5.26 (original: Q26)
std::array<q<0,31>, 256> speakerBuffer; // Q0.31 (original: Q31)
std::array<q<5,26>, 256> irEstimate;    // Q5.26 (original: Q26)
float muFinal;

for (uint32_t idx = 0; idx < config.filterLength; idx++) {
    // Gradient estimate: Q5.26 * Q0.31 -> Q5.26 (explicit format)
    auto gradientEstimate = mul_as<5,26>(error, speakerBuffer[idx]);

    // Convert muFinal to fixed-point for the calculation
    q<0,31> muFinal_fixed = q<0,31>::from_float(muFinal);

    // Filter step: Q5.26 * Q0.31 -> Q5.26 (explicit format, fits in 32 bits)
    auto filterStep = mul_as<5,26>(gradientEstimate, muFinal_fixed);

    // Addition with automatic saturation (Q5.26 + Q5.26 -> Q5.26)
    irEstimate[idx] = irEstimate[idx] + filterStep;
    // Saturation happens automatically when result exceeds Q5.26 range!
}
```

**Benefits:**
- ✅ Automatic saturation on conversion
- ✅ No manual min/max clamping
- ✅ Type system prevents accidental format mismatches
- ✅ Clear conversion from float to fixed-point
- ✅ Explicit output formats for multiplications

---

## Summary of Advantages

### Type Safety
- **Before:** Q-formats tracked in comments, easy to get wrong
- **After:** Q-formats are part of the type system, compiler-verified

### Error Prevention
- **Before:** Manual shifts, saturation, and rounding - prone to off-by-one errors
- **After:** Library handles all low-level details correctly

### Readability
- **Before:** Magic numbers (31, 26, 5) scattered throughout code
- **After:** Self-documenting types like `q<5,26>` and `mul_as<5,26>()`

### Maintainability
- **Before:** Changing Q-format requires updating shifts throughout the codebase
- **After:** Change the type declaration, compiler ensures consistency

### Performance
- **Before:** Manual implementation may miss optimization opportunities
- **After:** Backend system (Reference, Xtensa) provides optimized implementations automatically

### Portability
- **Before:** HiFi3 intrinsics mixed with business logic
- **After:** Business logic is backend-agnostic, backends handle optimization

---

## Scalar vs Vectorized Operations

The library provides both scalar and vectorized operations. Choose based on your use case:

### Scalar Operations (Single Values)
```cpp
// Use * operator or mul_as<> for single sample operations
q<5,26> a, b;
auto result = mul_as<5,26>(a, b);  // Single multiplication
```

**When to use**: Sample-by-sample processing (e.g., energy calculation per sample)

### Vectorized Operations (Arrays)
```cpp
// Use vec_elemult for element-wise array operations
q<5,26> array_a[256], array_b[256], result[256];
fp::vec_elemult<5,26>(result, array_a, array_b, 256);  // Vectorized

// Use specialized functions for reductions
auto dot_result = fp::dot<5,26>(array_a, array_b, 256);  // Dot product
auto mean_result = fp::mean(array_a, 256);                // Mean
```

**When to use**: Batch processing of arrays

**Why it matters**: Vectorized operations leverage SIMD instructions (e.g., NatureDSP on Xtensa) for significant performance improvements.

**Important Distinction**:

| Operation | Precision | Rounding | Implementation | Use Case |
|-----------|-----------|----------|----------------|----------|
| `mimi_mul_32x32` | 64-bit intermediate | None (truncate) | Software scalar | Single sample |
| `vec_elemult32x32` (NatureDSP) | 32-bit only | Yes | Hardware SIMD | Array processing |
| `mul_as<>()` (fp_lib) | 64-bit intermediate | Yes (round-to-nearest) | Software scalar | Single sample |
| `fp::vec_elemult<>()` (fp_lib + Xtensa) | 32-bit only | Yes (backend-specific) | Hardware SIMD | Array processing |

**Code Example**:
```cpp
// Original mimi_mul_32x32 - scalar with 64-bit intermediate, NO ROUNDING
inline int32_t mimi_mul_32x32(int32_t x, int32_t y, uint32_t shift) {
    // trust the caller that the result won't overflow
    return (int32_t) (((int64_t) x * y) >> shift);  // Just truncate!
}

// fp_lib scalar equivalent (sample-by-sample)
q<5,26> a, b;
auto result = mul_as<5,26>(a, b);  // Uses 64-bit intermediate + ROUNDING

// fp_lib vectorized (array processing, uses NatureDSP on Xtensa)
q<5,26> arr_a[256], arr_b[256], result[256];
fp::vec_elemult<5,26>(result, arr_a, arr_b, 256);  // 32-bit SIMD with rounding!
```

**Important Note on Rounding**:
- **`mimi_mul_32x32`**: No rounding - simple right shift (truncation)
- **`mul_as<>()`**: Uses round-to-nearest (via `round_shift` in helpers.hpp)
- **`fp::vec_elemult<>()`**: Backend-dependent rounding (NatureDSP uses specific rounding mode)

This means results may differ slightly between manual implementation and fp_lib. The fp_lib approach provides better numerical accuracy with rounding.

**Key Takeaway**: For array operations, use `fp::vec_elemult<>()` to get the NatureDSP 32-bit SIMD acceleration on Xtensa, which is much faster than scalar operations.

---

## Migration Strategy

1. **Identify Q-formats:** Document the Q-format of each variable (may require analysis)
2. **Replace buffer types:** Change `int32_t*` to `q<I,F>*`
3. **Replace operations:**
   - Change `mimi_mul_32x32()` to `*` or `mul_as<>()`
   - Use built-in functions: `fp::dot()`, `fp::mean()`, `fp::sum()`, etc.
4. **Remove manual shifts:** Use `.convert<I,F>()` for explicit conversions
5. **Remove saturation code:** Library handles it automatically
6. **Test incrementally:** Verify each function maintains bit-exact behavior
7. **Optimize:** Leverage backend system for platform-specific optimizations

## Next Steps

To complete this refactoring:
1. Create a wrapper that adapts the original function signatures
2. Implement helper types for common Q-formats used in NLMS
3. Add unit tests comparing old vs new implementation
4. Profile to ensure performance is maintained or improved
5. Gradually migrate related functions (updateEnvelope, etc.)
   - Note: `sum` and `mean` are already provided by the library!
