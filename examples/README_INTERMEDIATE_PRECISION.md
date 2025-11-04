# Intermediate Precision Example

## Overview

This example demonstrates the critical differences between scalar and vectorized fixed-point operations, with a focus on **intermediate precision** requirements for Q-format arithmetic.

## The Problem

When converting code like `mimi_mul_32x32` to vectorized operations, you might assume you can simply use NatureDSP's `vec_elemult32x32`. **This is wrong and will produce incorrect results!**

## Why?

### Scalar Operation (mimi_mul_32x32)
```cpp
int32_t mimi_mul_32x32(int32_t x, int32_t y, uint32_t shift) {
    return (int32_t)(((int64_t) x * y) >> shift);  // 64-bit intermediate!
}
```
- Uses **64-bit intermediate** (int64_t)
- Can hold full 32×32 product without overflow
- Works correctly for Q-format arithmetic

### Vectorized Element-wise (vec_elemult32x32)
```cpp
void vec_elemult32x32(int32_t *z, int32_t *x, int32_t *y, int N);
// Performs: z[i] = SAT32(x[i] * y[i])  // Only 32-bit!
```
- Uses **32-bit intermediate** only
- Product overflows and saturates
- **Incorrect for Q-format arithmetic**

### Example: Q0.31 × Q0.31

```
Value: 0.5 × 0.5 = 0.25
Q0.31: 1073741824 × 1073741824 = 1152921504606846976

With 64-bit intermediate:
  1152921504606846976 >> 31 = 536870912 (0.25 in Q0.31) ✓

With 32-bit intermediate:
  1073741824 × 1073741824 = OVERFLOW → INT32_MAX
  INT32_MAX is NOT 0.25 ✗
```

## When Can You Vectorize?

| Operation | NatureDSP Function | Intermediate | Works for Q-format? |
|-----------|-------------------|--------------|---------------------|
| Element-wise multiply | `vec_elemult32x32` | 32-bit | **NO** ✗ |
| Dot product | `vec_dot32x32` | 64-bit accum | **YES** ✓ |
| Sum | `vec_sum32x32` | 64-bit accum | **YES** ✓ |
| Power (sum of squares) | `vec_power32x32` | 64-bit accum | **YES** ✓ |

## fp_lib Solution

The library **automatically** chooses the correct implementation:

```cpp
// Element-wise multiply: Falls back to software (64-bit intermediate)
fp::array_elemult<5,26>(output, a, b, count);
// Uses reference implementation, not vec_elemult32x32

// Dot product: Uses hardware SIMD (64-bit accumulator)
auto result = fp::dot<5,26>(a, b, count);
// Uses vec_dot32x32 on Xtensa for performance
```

## Refactoring Guide

### 1. Sample-by-Sample Multiplication
```cpp
// Before:
energy[i] = mimi_mul_32x32(sample, sample, 36);

// After:
energy[i] = mul_as<5,26>(sample, sample);
```

### 2. Array Element-wise Multiplication
```cpp
// Before:
for (i = 0; i < N; i++) {
    output[i] = mimi_mul_32x32(a[i], b[i], shift);
}

// After (still uses software, but type-safe):
fp::array_elemult<5,26>(output, a, b, N);
// Note: Can't vectorize this on Xtensa!
```

### 3. Dot Product (Can Vectorize!)
```cpp
// Before:
int32_t acc = 0;
for (i = 0; i < N; i++) {
    acc += mimi_mul_32x32(a[i], b[i], shift);
}

// After (uses SIMD!):
auto acc = fp::dot<5,26>(a, b, N);
// Uses vec_dot32x32 on Xtensa
```

### 4. Sum of Squares (Can Vectorize!)
```cpp
// Before:
int32_t power = 0;
for (i = 0; i < N; i++) {
    power += mimi_mul_32x32(a[i], a[i], shift);
}

// After (uses SIMD!):
auto power = fp::power<5,26>(a, N);
// Uses vec_power32x32 on Xtensa
```

## Running the Example

```bash
cd examples
g++ -std=c++17 -I.. intermediate_precision_example.cpp -o intermediate_precision
./intermediate_precision
```

Or with CMake:
```bash
cmake --build build --target intermediate_precision_example
./build/intermediate_precision_example
```

## Key Takeaways

1. **mimi_mul_32x32 cannot be directly vectorized** with `vec_elemult32x32`
2. **Element-wise multiplication** requires 64-bit intermediate for Q-format
3. **Reduction operations** (dot, sum, power) CAN be vectorized with NatureDSP
4. **fp_lib automatically** uses the right implementation
5. **Performance**: Vectorized reductions are ~10-20x faster than scalar loops

## See Also

- `porting/fixedpoint_refactoring_example.md` - General refactoring guide
- `backends/xtensa/vector_elemwise.hpp` - Implementation details
- `backends/xtensa/vector_ops.hpp` - Vectorized reduction operations
- `README-NatureDSP.md` - NatureDSP function reference
