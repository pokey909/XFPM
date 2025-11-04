# Fixed-Point Library Examples

## Overview

This directory contains educational examples demonstrating the fixed-point library:

1. **basic_tutorial.cpp** - Introduction to scalar and vector operations
2. **vectorized_operations.cpp** - Advanced vectorization and NatureDSP optimization
3. **intermediate_precision_example.cpp** - Understanding intermediate precision in scalar vs vectorized operations

---

## Basic Tutorial (`basic_tutorial.cpp`)

Comprehensive educational example demonstrating:

### Part 1: Basic Operations - Same Precision
- Multiplication and division in Q1.15 format
- Saturation handling (division result overflow demonstration)
- Using `div_as<>` for explicit output format control

### Part 2: Mixed Precision - Transparent Requantization
- Mixing Q4.12 and Q1.15 formats seamlessly
- Automatic requantization by the library
- Explicit format control with `mul_as<>` and `div_as<>`

**Key Takeaway:** No manual bit shifting or requantization needed!

### Part 3: Scalar Operations (Element-wise)
Demonstrates single-value operations:
- **Logarithms:** `log2()`, `logn()`, `log10()`
- **Square root:** `sqrt()`
- **Power:** `pow(exponent)`
- **Trigonometric:** `sin()`, `cos()`, `tan()`, `atan()`
- **Activation:** `sigmoid()`

### Part 4: Vector Operations (Array-based)
Array-level operations with two API styles:

**Option 1 - FixedPointArray wrapper:**
```cpp
FixedPointArray<1, 15> arr(data, length);
auto min_val = arr.min();
auto max_val = arr.max();
auto sum_val = arr.sum();
```

**Option 2 - Static methods:**
```cpp
auto min_val = Q1_15::array_min(data, length);
auto max_val = Q1_15::array_max(data, length);
```

**In-place operations:**
- `arr.scale(factor)` - Multiply all elements by a scalar
- `arr.shift(amount)` - Bit-shift all elements

**Vector operations:**
- `arr.dot_product(other)` - Dot product of two arrays (Note: current implementation has precision limitations)

### Part 5: Complex Example - Neural Network Layer
Shows composition of operations:
1. Dot product for weighted sum
2. Sigmoid activation function
3. Demonstrates building blocks for ML inference

### Part 6: Backend Selection
- **ReferenceBackend:** Portable C++ implementation
- **XtensaBackend:** Optimized for Cadence Xtensa HiFi DSP

Both backends provide identical API - only performance differs!

---

## Vectorized Operations Tutorial (`vectorized_operations.cpp`)

Advanced tutorial focusing on array operations and hardware optimization:

### Part 1: Scalar vs Vectorized Operations
- Comparison of element-wise vs SIMD-style processing
- Demonstrates performance benefits of vectorization
- Shows when to use scalar vs vector operations

### Part 2: NatureDSP Fast Variants - Restrictions
Comprehensive documentation of optimization variants:

**Operations with fast variants:**
- ✓ `vec_dot16x16_fast` / `vec_dot32x32_fast` - Dot product
- ✓ `vec_add16x16_fast` / `vec_add32x32_fast` - Element-wise addition
- ✓ `vec_power16x16_fast` / `vec_power32x32_fast` - Sum of squares

**Fast variant requirements (ALL must be met):**
1. All pointers must be **8-byte aligned**
2. Array length N must be **≥ 4**
3. Array length N must be a **multiple of 4**

**Automatic fallback:** If any requirement is not satisfied, the library automatically falls back to regular variants.

### Part 3: Alignment Requirements Demonstration
Three detailed scenarios:

**Scenario 1: Properly Aligned** ✓
- All pointers 8-byte aligned
- Length = 8 (multiple of 4)
- Result: Uses fast variant

**Scenario 2: Misaligned Pointer** ✗
- One pointer not aligned
- Falls back to regular variant
- Still works correctly!

**Scenario 3: Wrong Length** ✗
- All pointers aligned
- Length = 7 (not multiple of 4)
- Falls back to regular variant

### Part 4: Vectorized Operations - Complete Examples

**Element-wise operations:**
```cpp
arr_a.add(arr_b, arr_out);      // vec_add16x16_fast
arr_a.sub(arr_b, arr_diff);     // vec_elesub16x16
arr_a.elemult(arr_b, arr_prod); // vec_elemult + vec_shift
```

**Reduction operations:**
```cpp
auto sum = arr.sum();           // vec_sum16x16
auto min = arr.min();           // vec_min16x16
auto max = arr.max();           // vec_max16x16
```

**Dot product:**
```cpp
auto dot = v1.dot_product(v2);  // vec_dot16x16_fast
```

**Statistical operations:**
```cpp
auto mean = arr.mean();         // vec_mean16x16
auto rms = arr.rms();           // vec_rms16x16
auto variance = arr.variance(); // vec_var16x16
auto stddev = arr.stddev();     // vec_stddev16x16
```

### Part 5: Performance Considerations

**When to use vectorized operations:**
- ✓ Processing large arrays (100s-1000s of elements)
- ✓ Same operation on all elements
- ✓ Data can be aligned to 8 bytes
- ✓ Known lengths that are multiples of 4

**When NOT to use:**
- ✗ Only 1-2 elements (scalar is simpler)
- ✗ Each element needs different operations
- ✗ Inherently misaligned data

**Tips for maximum performance:**

1. **Alignment:**
   ```cpp
   alignas(8) int16_t my_array[N];
   ```

2. **Dynamic allocation:**
   ```cpp
   void* ptr = aligned_alloc(8, N * sizeof(int16_t));
   ```

3. **Use lengths that are multiples of 4:** 4, 8, 12, 16, 20, ...

**Expected speedup on Xtensa HiFi4:**
- Dot product: 2-3× faster
- Element-wise add: 2-4× faster
- Power (sum of squares): 1.5-2× faster

### Part 6: Backend Comparison
Side-by-side comparison of Reference vs Xtensa backends:
- Identical API
- Identical results
- Different performance characteristics

### Part 7: Real-World Example - Audio Processing
Complete audio mixer implementation:
1. Apply channel gains (element-wise multiply)
2. Mix channels (element-wise add)
3. Compute RMS level (for VU meter)

**All operations use fast variants** when data is properly aligned!

---

## Intermediate Precision Example (`intermediate_precision_example.cpp`)

**Critical guide for understanding when vectorization works for Q-format arithmetic!**

This example answers the question: *"Can I replace `mimi_mul_32x32` with NatureDSP's `vec_elemult32x32` for performance?"*

**Short answer: NO for element-wise multiplication, YES for reductions (dot product, sum)!**

### Part 1: Scalar Operations - 64-bit Intermediate

Demonstrates `mimi_mul_32x32` style operations:
- Uses **64-bit intermediate** (int64_t) for products
- Prevents overflow in Q-format arithmetic
- Example: Q0.31 × Q0.31 needs 64-bit to hold result before shift
- fp_lib equivalent: `mul_as<>()` (adds rounding for better accuracy)

**Key insight:** Manual scalar multiplication is correct because it uses 64-bit intermediate!

### Part 2: Vectorized Element-wise - NO 64-bit Intermediate ⚠️

**WARNING:** NatureDSP `vec_elemult32x32` is **NOT** suitable for Q-format!

```cpp
void vec_elemult32x32(int32_t *z, int32_t *x, int32_t *y, int N);
// Performs: z[i] = SAT32(x[i] * y[i])  ← Only 32-bit intermediate!
```

**Why it fails:**
- Q0.31 × Q0.31 produces Q0.62 intermediate (needs 64 bits)
- `vec_elemult32x32` uses only 32-bit intermediate
- Result: **Saturation to INT32_MAX** instead of correct value
- Example: 0.5 × 0.5 = 0.25 becomes SAT (wrong!)

**fp_lib solution:**
- Automatically detects this limitation
- Falls back to reference implementation (software loop with 64-bit intermediate)
- Correctness > Performance for element-wise multiplication

### Part 3: Vectorized Reduction - 64-bit Accumulator ✓

**GOOD NEWS:** NatureDSP reduction operations CAN be used!

```cpp
int64_t vec_dot32x32(const int32_t *x, const int32_t *y, int N);
// Returns: int64_t accumulator ← Full 64-bit precision!
```

**Why it works:**
- Each product: x[i] × y[i] produces 64-bit intermediate
- Accumulator: 64-bit register holds sum
- User applies final shift to convert Q-format
- Result: **Correct and FAST!**

**Vectorizable operations:**
- ✓ Dot product (`vec_dot32x32`)
- ✓ Sum (`vec_sum32x32`)
- ✓ Power/Energy (`vec_power32x32` - sum of squares)
- ✓ Mean (uses `vec_sum32x32` internally)

### Part 4: Comparison Table

Comprehensive table showing:
- mimi_mul_32x32: 64-bit intermediate (truncate) → Use `mul_as<>()`
- vec_elemult32x32: 32-bit only (saturates) → **Cannot use for Q-format!**
- vec_dot32x32: 64-bit accumulator → **Perfect for Q-format!**
- vec_sum32x32: 64-bit accumulator → **Perfect for Q-format!**

### Part 5: Practical Guidance

**When refactoring from `mimi_mul_32x32`:**

1. **Sample-by-sample multiplication:**
   ```cpp
   // Before:
   energy[i] = mimi_mul_32x32(sample, sample, 36);

   // After:
   energy[i] = fp::mul_as<5,26>(sample, sample);  // Still scalar
   ```

2. **Element-wise array multiplication:**
   ```cpp
   // Before:
   for (i = 0; i < N; i++) {
       output[i] = mimi_mul_32x32(a[i], b[i], shift);
   }

   // After:  Still uses scalar loop (no SIMD for Q-format!)
   for (i = 0; i < N; i++) {
       output[i] = fp::mul_as<5,26>(a[i], b[i]);
   }
   ```

3. **Dot product (CAN vectorize!):**
   ```cpp
   // Before:
   for (i = 0; i < N; i++) {
       acc += mimi_mul_32x32(a[i], b[i], shift);
   }

   // After:  Uses vec_dot32x32 on Xtensa!
   auto acc = fp::dot<5,26>(a, b, N);  // 10-20× faster!
   ```

4. **Sum of squares (CAN vectorize!):**
   ```cpp
   // Before:
   for (i = 0; i < N; i++) {
       power += mimi_mul_32x32(a[i], a[i], shift);
   }

   // After:  Uses vec_power32x32 on Xtensa!
   auto power = fp::power<5,26>(a, N);  // Much faster!
   ```

### Part 6: Key Takeaways

1. **mimi_mul_32x32 cannot be directly vectorized** with `vec_elemult32x32`
2. **Element-wise multiplication** requires 64-bit intermediate → Must use scalar loop
3. **Reduction operations** (dot, sum, power) CAN vectorize → Huge performance gain
4. **fp_lib automatically** chooses correct implementation
5. **Performance**: Vectorized reductions are 10-20× faster than scalar loops

### Part 7: Real-World Scenarios

Shows which NLMS algorithm operations can/cannot vectorize:
- Energy per sample: **Cannot vectorize** (element-wise multiply)
- Dot product (FIR): **CAN vectorize** (reduction)
- Sum of energy: **CAN vectorize** (reduction)
- Mean calculation: **CAN vectorize** (uses sum)

## Building and Running

```bash
# Build all examples
cmake --build build

# Or build individually
cmake --build build --target basic_tutorial
cmake --build build --target vectorized_operations
cmake --build build --target intermediate_precision_example

# Run the basic tutorial
./build/basic_tutorial

# Run the vectorized operations tutorial
./build/vectorized_operations

# Run the intermediate precision example
./build/intermediate_precision_example
```

## Key Concepts Demonstrated

### 1. Type Safety
All operations are type-checked at compile time. Mixing incompatible formats produces compiler errors.

### 2. Automatic Precision Management
```cpp
using Q4_12 = q<4, 12>;  // [-8, 8) range
using Q1_15 = q<1, 15>;  // [-1, 1) range

auto a = Q4_12::from_float(3.5f);
auto b = Q1_15::from_float(0.25f);
auto result = a * b;  // Library handles all shifting automatically!
```

### 3. Explicit vs. Implicit Format Control
```cpp
// Implicit: result takes format of left operand
auto result1 = a * b;  // Q4.12 format

// Explicit: specify output format
auto result2 = mul_as<3, 13>(a, b);  // Q3.13 format
```

### 4. Saturation Awareness
When results exceed the output format's range, they saturate:
```cpp
using Q1_15 = q<1, 15>;  // Can only hold [-1, 1)
auto a = Q1_15::from_float(0.75f);
auto b = Q1_15::from_float(0.5f);
auto result = a / b;  // Would be 1.5, saturates to ~1.0

// Solution: use larger output format
auto correct = div_as<3, 13>(a, b);  // Now holds 1.5 correctly
```

### 5. Backend Abstraction
```cpp
// Change backend by changing template parameter
using Q16_ref = q<1, 15, ReferenceBackend>;
using Q16_xtensa = q<1, 15, XtensaBackend>;

// API is identical for both!
auto result1 = a_ref.sqrt();
auto result2 = a_xtensa.sqrt();
```

## Common Patterns

### Signal Processing
```cpp
auto signal = Q1_15::from_float(audio_sample);
auto filtered = signal * filter_coef;
auto output = filtered.sqrt();  // RMS calculation
```

### Control Systems
```cpp
auto error = setpoint - measured_value;
auto control_signal = error * kp;  // Proportional control
```

### Neural Networks
```cpp
auto activation = weighted_sum.sigmoid();
auto relu_out = value.relu();
```

## Notes

- **Precision Loss:** Some vector operations (like `dot_product`) use raw integer multiplication and may have precision limitations in the current implementation
- **Saturation:** Always ensure your Q format can hold the expected range of results
- **Performance:** XtensaBackend provides hardware-accelerated operations when running on Xtensa DSP cores

## See Also

- `fp.hpp` - Main library header with full API documentation
- `tests/` - Comprehensive test suite with more examples
- `backends/` - Backend implementations (Reference and Xtensa)
