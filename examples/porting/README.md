# Fixed-Point Library Porting Examples

This directory contains real-world examples demonstrating how to refactor manual fixed-point arithmetic code to use the `fp_lib` type-safe library.

## Overview

The examples are based on the NLMS (Normalized Least Mean Squares) adaptive filter algorithm from `nlms_node.cpp`, which contains extensive manual fixed-point arithmetic with:
- Q-format tracking via comments (Q26, Q31, etc.)
- Manual bit shifts for format conversions
- Manual saturation logic
- Helper functions like `mimi_mul_32x32`
- Platform-specific intrinsics (HiFi3)

These examples show how `fp_lib` eliminates these pain points while improving type safety, readability, and maintainability.

## Files

### 1. [fixedpoint_refactoring_example.md](fixedpoint_refactoring_example.md)
**Comprehensive refactoring guide with multiple examples**

Contains five detailed examples:
- **Dot Product**: Shows basic operation using built-in `fp::dot()`
- **Energy Calculation**: Demonstrates format conversion handling
- **Complex Processing**: Shows multi-step operations
- **Mean Calculation**: Demonstrates built-in `fp::mean()` function
- **Gradient Update**: Demonstrates saturation handling

Each example includes:
- Before/after code comparison
- Issues with manual approach
- Benefits of fp_lib approach
- Migration strategy

📖 **Start here** for a complete overview of the refactoring process.

### 2. [dot_product_comparison.cpp](dot_product_comparison.cpp)
**Compilable side-by-side comparison**

Features:
- Dot product implementation (both versions)
- Energy calculation example
- Full working code that compiles and runs
- Demonstrates type safety benefits

```bash
# To compile and run:
g++ -std=c++17 -I../.. dot_product_comparison.cpp -o dot_comparison
./dot_comparison
```

💡 **Best for**: Seeing concrete, runnable examples

### 3. [process_sample_refactored.cpp](process_sample_refactored.cpp)
**Complete function refactoring**

Shows full transformation of `AfcNlmsNode::processSample`:
- Complete class with all buffers and state
- Before: 80+ lines with manual fixed-point
- After: ~60 lines with fp_lib, much clearer
- Demonstrates automatic saturation
- Shows type aliases for complex formats

```bash
# To compile and run:
g++ -std=c++17 -I../.. process_sample_refactored.cpp -o process_sample
./process_sample
```

🎯 **Best for**: Understanding large-scale refactoring

### 4. [nlms_node.cpp](nlms_node.cpp)
**Original source code (reference only)**

The original NLMS implementation showing:
- Float version (lines 144-202)
- HiFi3 optimized fixed-point (lines 247-318)
- Reference fixed-point (lines 351-432)

⚠️ **Reference material**: Study this to understand the problem being solved

## Key Advantages of fp_lib

### Type Safety
```cpp
// Before: Q-format only in comments
int32_t buffer[256];  // Q26 - hope we remember!

// After: Q-format in type system
std::array<q<5, 26>, 256> buffer;  // Compiler enforced!
```

### No Magic Numbers
```cpp
// Before: Manual shift calculation
result = mimi_mul_32x32(a, b, 31 + 5);  // Why 31 + 5?

// After: Explicit format
result = mul_as<10, 52>(a, b);  // Clear output format!
```

### Automatic Saturation
```cpp
// Before: 9 lines of manual saturation
static constexpr int32_t q26Max = (1 << 26) - 1;
static constexpr int32_t q26Min = -(1 << 26);
int32_t irSample = ir_estimate[idx] + filterStep;
if (irSample < q26Min) {
    irSample = q26Min;
} else if (irSample > q26Max) {
    irSample = q26Max;
}
ir_estimate[idx] = irSample;

// After: 1 line with automatic saturation
ir_estimate[idx] = ir_estimate[idx] + filterStep.convert<1, 30>();
```

### Self-Documenting
```cpp
// Before: Need comments to track format
int32_t energy;  // Q26 format

// After: Type documents itself
q<10, 52> energy;  // 10 integer bits, 52 fractional bits
```

### Easy Refactoring
To change Q-format from Q5.26 to Q8.23:

**Before**: Find and update dozens of shift values throughout code
```cpp
result = mimi_mul_32x32(a, b, 31 + 5);  // Need to recalculate
result = result >> 5;                   // Need to update
```

**After**: Change one type declaration
```cpp
using SpeakerFormat = q<8, 23>;  // Done! Everything else updates automatically
```

### Backend Portability
```cpp
// Same code works with different backends
template<typename Backend = ReferenceBackend>
class AfcNlmsNode {
    FixedPoint<5, 26, Backend> data;  // Uses ReferenceBackend
};

// Switch to Xtensa backend for HiFi3 optimization:
AfcNlmsNode<XtensaBackend> node;  // Automatic SIMD optimization!
```

## Migration Workflow

1. **Analyze**: Document Q-formats used in existing code
2. **Replace Types**: Change `int32_t*` to `q<I,F>*`
3. **Replace Operations**:
   - `mimi_mul_32x32(a, b, shift)` → `a * b` or `mul_as<I,F>(a, b)`
   - Manual shifts → `.convert<I,F>()`
4. **Remove Saturation**: Library handles it automatically
5. **Test**: Verify bit-exact behavior
6. **Optimize**: Leverage backend system

## Quick Start

1. Read [fixedpoint_refactoring_example.md](fixedpoint_refactoring_example.md) for concepts
2. Compile and run [dot_product_comparison.cpp](dot_product_comparison.cpp)
3. Study [process_sample_refactored.cpp](process_sample_refactored.cpp) for complex example
4. Apply to your own code!

## Performance Notes

- **Reference Backend**: Pure C++, portable, good for development
- **Xtensa Backend**: HiFi3-optimized, leverages NatureDSP library
- Switching backends requires only changing template parameter
- No performance penalty vs. manual implementation (often faster)

## Related Documentation

- [../../README-NatureDSP.md](../../README-NatureDSP.md) - NatureDSP function reference
- [../../CLAUDE.md](../../CLAUDE.md) - Project architecture guide
- [../../fp.hpp](../../fp.hpp) - Core library documentation

## Built-in Library Functions

The library provides many operations you don't need to implement yourself:

### Vector Operations
```cpp
// Dot product
auto result = fp::dot(buffer_a, buffer_b, count);

// Sum
auto total = fp::sum(buffer, count);

// Mean (average)
auto average = fp::mean(buffer, count);

// Element-wise multiplication (vectorized, replaces mimi_mul_32x32 for arrays)
fp::vec_elemult<OutI, OutF>(result, buffer_a, buffer_b, count);

// Power (raise all elements)
fp::power(buffer, count, exponent);

// Batch operations with explicit output format
auto result = fp::dot<OutI, OutF>(buffer_a, buffer_b, count);
```

### Mathematical Operations
```cpp
// Square root
auto root = fp::sqrt(value);

// Logarithm/antilogarithm
auto log_val = fp::log2(value);
auto exp_val = fp::antilog2(value);

// Division
auto quotient = fp::divide(numerator, denominator);
```

### Activation Functions
```cpp
// Sigmoid, ReLU, softmax, tanh
auto activated = fp::sigmoid(value);
```

All of these operations:
- ✅ Are backend-optimized (SIMD when available)
- ✅ Handle Q-format conversions automatically
- ✅ Provide type-safe interfaces
- ✅ Support explicit output format specification

**Important Note on Scalar vs Vectorized Operations**:
- **`mimi_mul_32x32`**: Scalar multiplication with **64-bit intermediate** precision, **NO rounding** (truncation only)
  ```cpp
  // Original - just truncates!
  return (int32_t) (((int64_t) x * y) >> shift);
  ```
- **`vec_elemult32x32` (NatureDSP)**: Vectorized SIMD working entirely in **32-bit** with hardware acceleration and rounding
- **For refactoring**:
  - Single sample operations: Use `mul_as<>()` (64-bit intermediate, **with rounding** - slightly different from `mimi_mul_32x32`)
  - Array operations: Use `fp::vec_elemult()` (32-bit SIMD with Xtensa backend, like NatureDSP, with rounding)

**Rounding Behavior**:
- Original `mimi_mul_32x32`: Truncates (no rounding)
- fp_lib `mul_as<>()`: Round-to-nearest via `round_shift`
- Results may differ slightly, but fp_lib provides better numerical accuracy

## Questions?

Common patterns:
- **Built-in operations**: `fp::dot()`, `fp::mean()`, `fp::sum()`, `fp::sqrt()`, etc.
- **Multiplication**: `mul_as<OutI, OutF>(a, b)` for explicit format
- **Format conversion**: `value.convert<NewI, NewF>()`
- **From float**: `q<I, F>::from_float(f)`
- **To float**: `value.to_float()`
- **Type aliases**: `using MyFormat = q<I, F, Backend>;`

See the examples for complete patterns!
