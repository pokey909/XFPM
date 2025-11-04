# Xtensa HiFi3 Cross-Compilation Guide

This project supports building both native (x86-64) and Xtensa HiFi3 executables from a single CMake configuration.

## Configuration

**Xtensa Toolchain:**
- Path: `${HOME}/xtensa/XtDevTools/install/tools/RJ-2023.2-linux`
- Core: `HiFi3_Aria2_0_RJ2023_2`
- Compiler: `xt-clang++`

## Building

### Native Build Only (Default)
```bash
cmake -B build -S .
cmake --build build
```

### Native + Xtensa Build
```bash
cmake -B build -S . -DBUILD_XTENSA=ON
cmake --build build
```

This creates:
- **7 native executables** (x86-64): `test_multiply`, `test_divide`, etc.
- **7 Xtensa executables** (Tensilica Xtensa): `test_multiply_xtensa`, `test_divide_xtensa`, etc.

## Running Tests

### Native Tests Only
```bash
cd build && ctest -R "^(Multiply|Divide|Logarithm|Antilogarithm|SquareRoot|Power|ArrayOperations)$"
```

### All Tests (14 total)
```bash
cd build && ctest -N  # List all tests
```

**Note:** Xtensa tests (`*_Xtensa`) cannot run on x86-64 host. They must be transferred to Xtensa hardware or executed in an ISS (Instruction Set Simulator).

## Binary Verification

Verify executables are compiled for correct architecture:
```bash
file build/test_multiply          # Should show: x86-64
file build/test_multiply_xtensa   # Should show: Tensilica Xtensa
```

## Test List

### Native Tests (7)
1. Multiply
2. Divide
3. Logarithm
4. Antilogarithm
5. SquareRoot
6. Power
7. ArrayOperations

### Xtensa Tests (7)
8. Multiply_Xtensa
9. Divide_Xtensa
10. Logarithm_Xtensa
11. Antilogarithm_Xtensa
12. SquareRoot_Xtensa
13. Power_Xtensa
14. ArrayOperations_Xtensa

## Toolchain File

A standalone Xtensa toolchain file is available at:
```
cmake/xtensa-toolchain.cmake
```

This can be used for full cross-compilation projects:
```bash
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=cmake/xtensa-toolchain.cmake
```
