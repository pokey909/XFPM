# Xtensa HiFi3 Toolchain File
# Based on Cadence Xtensa Tools RJ-2023.2

set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR xtensa)

# Xtensa toolchain paths
set(XTENSA_TOOLS_ROOT "/home/alenhardt/xtensa/XtDevTools/install/tools/RJ-2023.2-linux")
set(XTENSA_CORE "HiFi3_Aria2_0_RJ2023_2")
set(XTENSA_SYSTEM "${XTENSA_TOOLS_ROOT}/XtensaTools/config")

# Xtensa build root for core-specific configurations
set(XTENSA_BUILD "${XTENSA_TOOLS_ROOT}/XtensaTools/Tools/build")

# Set compilers
set(CMAKE_C_COMPILER "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-clang")
set(CMAKE_CXX_COMPILER "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-clang++")
set(CMAKE_AR "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-ar")
set(CMAKE_RANLIB "${XTENSA_TOOLS_ROOT}/XtensaTools/bin/xt-ranlib")

# Core-specific flags
set(XTENSA_CORE_FLAGS "--xtensa-core=${XTENSA_CORE}")

# Compiler flags
set(CMAKE_C_FLAGS_INIT "${XTENSA_CORE_FLAGS} -mlongcalls")
set(CMAKE_CXX_FLAGS_INIT "${XTENSA_CORE_FLAGS} -mlongcalls")

# Don't search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Skip compiler checks (cross-compilation)
set(CMAKE_C_COMPILER_WORKS 1)
set(CMAKE_CXX_COMPILER_WORKS 1)

# Set Xtensa system config path
set(ENV{XTENSA_SYSTEM} "${XTENSA_SYSTEM}")
set(ENV{XTENSA_CORE} "${XTENSA_CORE}")
