CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(FixedPointLib)

# Option to build Xtensa executables
option(BUILD_XTENSA "Build Xtensa HiFi3 executables" OFF)

# Xtensa toolchain configuration
if(BUILD_XTENSA)
    set(XTENSA_TOOLS_ROOT "/home/alenhardt/xtensa/XtDevTools/install/tools/RJ-2023.2-linux")
    set(XTENSA_CORE "HiFi3_Aria2_0_RJ2023_2")
    set(XTENSA_CXX_COMPILER "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-clang++")
    set(XTENSA_CORE_FLAGS "--xtensa-core=${XTENSA_CORE}")

    # Set Xtensa environment variables
    set(ENV{XTENSA_SYSTEM} "${XTENSA_TOOLS_ROOT}/XtensaTools/config")
    set(ENV{XTENSA_CORE} "${XTENSA_CORE}")
endif()

# Function to add test executable with optional Xtensa version
function(add_test_executable target_name)
    # Parse arguments: source files
    set(source_files ${ARGN})

    # Add native executable
    add_executable(${target_name} ${source_files})
    target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    # Add Xtensa executable if enabled
    if(BUILD_XTENSA)
        set(xtensa_target "${target_name}_xtensa")
        set(xtensa_output "${CMAKE_BINARY_DIR}/${xtensa_target}")

        # Get absolute paths for source files
        set(abs_sources "")
        foreach(src ${source_files})
            get_filename_component(abs_src "${src}" ABSOLUTE)
            list(APPEND abs_sources "${abs_src}")
        endforeach()

        # Custom command to compile with Xtensa compiler
        add_custom_command(
            OUTPUT ${xtensa_output}
            COMMAND ${XTENSA_CXX_COMPILER}
                ${XTENSA_CORE_FLAGS}
                -mlongcalls
                -std=c++17
                -I${CMAKE_CURRENT_SOURCE_DIR}
                -I${CMAKE_CURRENT_SOURCE_DIR}/libs/ndsp/library/include
                ${abs_sources}
                ${CMAKE_CURRENT_SOURCE_DIR}/libs/ndsp/build/bin/NatureDSP_Signal-HiFi3_Aria2_0_RJ2023_2____mm0_release_clang-Xtensa-release.a
                -o ${xtensa_output}
            DEPENDS ${abs_sources}
            COMMENT "Building Xtensa executable: ${xtensa_target}"
            VERBATIM
        )

        # Custom target to trigger the build
        add_custom_target(${xtensa_target} ALL
            DEPENDS ${xtensa_output}
        )
    endif()
endfunction()

# Old monolithic test executable (kept for backwards compatibility)
add_executable(fp_demo main.cpp)

# Educational examples
add_executable(basic_tutorial examples/basic_tutorial.cpp)
target_include_directories(basic_tutorial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(vectorized_operations examples/vectorized_operations.cpp)
target_include_directories(vectorized_operations PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(intermediate_precision_example examples/intermediate_precision_example.cpp)
target_include_directories(intermediate_precision_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# AFC Float Observer Example
add_subdirectory(examples/afc/float)

# Debug test for vec_elemult behavior
add_test_executable(test_vec_elemult_direct test_vec_elemult_direct.cpp)

# New modular test suite (all tests in one executable)
add_executable(fp_tests
    tests/test_runner.cpp
    tests/test_multiply.cpp
    tests/test_divide.cpp
    tests/test_logarithm.cpp
    tests/test_antilogarithm.cpp
    tests/test_sqrt.cpp
    tests/test_power.cpp
    tests/test_array_ops.cpp
    tests/test_vector_ops.cpp
)
target_include_directories(fp_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(fp_tests PRIVATE FP_TEST_SUITE)

# Individual test executables
add_test_executable(test_multiply tests/test_multiply.cpp)
add_test_executable(test_divide tests/test_divide.cpp)
add_test_executable(test_logarithm tests/test_logarithm.cpp)
add_test_executable(test_antilogarithm tests/test_antilogarithm.cpp)
add_test_executable(test_sqrt tests/test_sqrt.cpp)
add_test_executable(test_power tests/test_power.cpp)
add_test_executable(test_array_ops tests/test_array_ops.cpp)
add_test_executable(test_vector_ops tests/test_vector_ops.cpp)

# Enable CTest support
enable_testing()

# Add individual tests to CTest (native)
add_test(NAME Multiply COMMAND test_multiply)
add_test(NAME Divide COMMAND test_divide)
add_test(NAME Logarithm COMMAND test_logarithm)
add_test(NAME Antilogarithm COMMAND test_antilogarithm)
add_test(NAME SquareRoot COMMAND test_sqrt)
add_test(NAME Power COMMAND test_power)
add_test(NAME ArrayOperations COMMAND test_array_ops)
add_test(NAME VectorOperations COMMAND test_vector_ops)

# Add Xtensa tests if enabled
if(BUILD_XTENSA)
    set(XT_RUN "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-run")
    add_test(NAME Multiply_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_multiply_xtensa)
    add_test(NAME Divide_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_divide_xtensa)
    add_test(NAME Logarithm_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_logarithm_xtensa)
    add_test(NAME Antilogarithm_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_antilogarithm_xtensa)
    add_test(NAME SquareRoot_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_sqrt_xtensa)
    add_test(NAME Power_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_power_xtensa)
    add_test(NAME ArrayOperations_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_array_ops_xtensa)
    add_test(NAME VectorOperations_Xtensa COMMAND ${XT_RUN} --xtensa-core=${XTENSA_CORE} ${CMAKE_BINARY_DIR}/test_vector_ops_xtensa)
endif()