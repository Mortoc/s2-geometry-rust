cmake_minimum_required(VERSION 3.18)
project(S2TestProgram)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(absl REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Threads REQUIRED)

# Add s2 as a subdirectory with tests disabled
set(BUILD_TESTS OFF)
set(BUILD_EXAMPLES OFF)
add_subdirectory(s2geometry-cpp)

# Create our test executable
add_executable(test_cpp_intersection test_cpp_intersection.cc)

# Link with s2 and required libraries
target_link_libraries(test_cpp_intersection s2)