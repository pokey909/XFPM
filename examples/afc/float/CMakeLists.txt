CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# AFC Float Observer Example
add_executable(afc_float_observer
    main.cpp
    nlms_node.cpp
)

# Include directories
target_include_directories(afc_float_observer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..  # For observer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../..  # For float/nlms_node.h
)

# Link math library if needed
target_link_libraries(afc_float_observer m)
